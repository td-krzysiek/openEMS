/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_cartesiansubgrid_conductivity.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include <stdexcept>

Engine_Ext_CartesianSubgrid_Conductivity::Engine_Ext_CartesianSubgrid_Conductivity(
	Operator_Ext_CartesianSubgrid_Conductivity *op)
	: Engine_Ext_CartesianSubgrid_Advance(op), m_Conductivity(op), m_Grid(op->m_Compact->Grid()),
	  m_OldElectric(m_Grid.EdgeCount()), m_SourceCurl(m_Grid.EdgeCount())
{
}

void Engine_Ext_CartesianSubgrid_Conductivity::Add(std::array<double, 2> &sum, double value)
{
	const double total = sum[0] + value;
	sum[1] += std::abs(sum[0]) >= std::abs(value) ? (sum[0] - total) + value : (value - total) + sum[0];
	sum[0] = total;
}

void Engine_Ext_CartesianSubgrid_Conductivity::Check(double timestep) const
{
	m_Grid.CheckStep(timestep);
	if (timestep != m_Conductivity->m_Timestep)
	{
		throw std::invalid_argument("CartesianSubgrid conductivity requires its factored timestep");
	}
}

void Engine_Ext_CartesianSubgrid_Conductivity::UpdateElectric(double timestep, double drive)
{
	Check(timestep);
	auto &g = m_Grid;
	const auto &c = *m_Conductivity;
	if (!c.m_HasLoss)
	{
		g.UpdateElectric(timestep, drive);
		return;
	}
	if (g.m_PendingMagneticTimestep != 0 || !std::isfinite(drive))
	{
		throw std::invalid_argument("CartesianSubgrid conductivity electric phase/drive is invalid");
	}
	std::copy(g.m_Electric.begin(), g.m_Electric.end(), m_OldElectric.begin());
	g.CurlTranspose(g.m_Magnetic, g.m_Increment);
	std::copy(g.m_Increment.begin(), g.m_Increment.end(), m_SourceCurl.begin());
	for (std::size_t i = 0; i < g.EdgeCount(); ++i)
	{
		g.m_Increment[i] -= c.m_Diagonal[i] * g.m_Electric[i];
	}
	for (const auto &term : c.m_OffDiagonal)
	{
		g.m_Increment[term.first.first] -= term.second * g.m_Electric[term.first.second];
		g.m_Increment[term.first.second] -= term.second * g.m_Electric[term.first.first];
	}
	c.Solve(g.m_Increment, m_Scratch);
	if (g.m_PortVoltage.empty())
	{
		for (std::size_t i = 0; i < g.EdgeCount(); ++i)
		{
			g.m_Electric[i] += timestep * g.m_Increment[i];
		}
	}
	else
	{
		g.ApplyResistors(timestep, drive, c.m_Responses, c.m_PortLower);
	}
	// Dissipation is evaluated on the implicit E midpoint BEFORE soft kicks.
	for (const auto &trace : c.m_Trace)
	{
		double field = 0;
		for (const auto &term : trace)
		{
			field += term.second * .5 * (m_OldElectric[term.first] + g.m_Electric[term.first]);
		}
		Add(m_ElectricLoss, timestep * field * field);
	}
	for (const auto &source : g.m_Sources)
	{
		const double change = drive * source.second;
		double derivative = g.m_Electric[source.first];
		for (const auto &term : g.m_MassOffDiagonal)
		{
			if (term.first.first == source.first)
			{
				derivative += term.second * g.m_Electric[term.first.second];
			}
			if (term.first.second == source.first)
			{
				derivative += term.second * g.m_Electric[term.first.first];
			}
		}
		g.m_SourceWork += change * (derivative + .5 * change - .5 * timestep * m_SourceCurl[source.first]);
		g.m_Electric[source.first] += change;
	}
	g.m_PendingMagneticTimestep = timestep;
}

void Engine_Ext_CartesianSubgrid_Conductivity::UpdateMagnetic(double timestep)
{
	Check(timestep);
	auto &g = m_Grid;
	const auto &c = *m_Conductivity;
	if (!c.m_HasLoss)
	{
		g.UpdateMagnetic(timestep);
		return;
	}
	if (g.m_PendingMagneticTimestep != timestep)
	{
		throw std::invalid_argument("CartesianSubgrid conductivity magnetic phase is invalid");
	}
	g.Curl(g.m_Electric, g.m_Curl);
	for (std::size_t i = 0; i < g.FaceCount(); ++i)
	{
		const double previous = g.m_Magnetic[i], rate = c.m_MagneticRate[i], a = .5 * timestep * rate;
		if (rate == 0)
		{
			g.m_Magnetic[i] -= timestep * g.m_Curl[i];
		}
		else
		{
			g.m_Magnetic[i] = ((1 - a) * previous - timestep * g.m_Curl[i]) / (1 + a);
			const double mid = .5 * (previous + g.m_Magnetic[i]);
			Add(m_MagneticLoss, timestep * rate * mid * mid);
		}
	}
	g.m_PendingMagneticTimestep = 0;
}

std::unique_ptr<GPU_CartesianSubgrid> Engine_Ext_CartesianSubgrid_Conductivity::CreateDevice(GPU_Backend &backend)
{
	if (!m_Conductivity->HasLoss())
	{
		return backend.CreateCartesianSubgrid(m_Grid);
	}
	const auto &c = *m_Conductivity;
	const auto index = [](std::size_t value)
	{
		if (value > std::numeric_limits<std::uint32_t>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: conductivity device index overflow");
		}
		return static_cast<std::uint32_t>(value);
	};
	const auto scalar = [](double value)
	{
		const float result = static_cast<float>(value);
		if (!std::isfinite(value) || !std::isfinite(result) ||
			(value != 0 && std::abs(result) < std::numeric_limits<float>::min()))
		{
			throw std::invalid_argument("CartesianSubgrid: unrepresentable conductivity device factor/coefficient");
		}
		return result;
	};
	auto data = m_Grid.DeviceData();
	if (data.rootConductance.empty())
	{
		data.portTimestep = 0;
	}
	data.blockStart = {0};
	data.blockEdge.clear();
	data.bandwidth.clear();
	data.lowerStart = {0};
	data.lower.clear();
	data.factorTails = {};
	const auto appendFactor = [&](const CartesianSubgrid::MassBlock &block)
	{
		for (std::size_t edge : block.edges)
		{
			data.blockEdge.push_back(index(edge));
		}
		data.blockStart.push_back(index(data.blockEdge.size()));
		data.bandwidth.push_back(index(block.halfBandwidth));
		const auto factor = CartesianSubgridFloat::Factor(block.lower, block.edges.size(),
			block.halfBandwidth, data.factorTails);
		data.lower.insert(data.lower.end(), factor.begin(), factor.end());
		data.lowerStart.push_back(index(data.lower.size()));
	};
	for (const auto &block : c.m_Factors)
	{
		appendFactor(block);
	}
	for (std::size_t edge = 0; edge < c.m_IndependentInverse.size(); ++edge)
	{
		if (c.m_IndependentInverse[edge] != 1)
		{
			CartesianSubgrid::MassBlock block;
			block.edges = {edge};
			block.halfBandwidth = 0;
			block.lower = {std::sqrt(1 / c.m_IndependentInverse[edge])};
			appendFactor(block);
		}
	}
	data.responseStart = {0};
	data.responseEdge.clear();
	data.responseWeight.clear();
	for (const auto &response : c.m_Responses)
	{
		for (const auto &term : response.weights)
		{
			data.responseEdge.push_back(index(term.first));
			data.responseWeight.push_back(scalar(term.second));
		}
		data.responseStart.push_back(index(data.responseEdge.size()));
	}
	data.portLower.clear();
	for (double value : c.m_PortLower)
	{
		data.portLower.push_back(scalar(value));
	}
	data.trackSourceWork = true;
	data.sourceOffStart = {0};
	data.sourceOffEdge.clear();
	data.sourceOffWeight.clear();
	for (const auto &source : m_Grid.m_Sources)
	{
		for (const auto &term : m_Grid.m_MassOffDiagonal)
		{
			if (term.first.first == source.first)
			{
				data.sourceOffEdge.push_back(index(term.first.second));
				data.sourceOffWeight.push_back(scalar(term.second));
			}
			if (term.first.second == source.first)
			{
				data.sourceOffEdge.push_back(index(term.first.first));
				data.sourceOffWeight.push_back(scalar(term.second));
			}
		}
		data.sourceOffStart.push_back(index(data.sourceOffEdge.size()));
	}
	CartesianSubgridConductivityData loss;
	loss.timestep = c.m_Timestep;
	loss.magneticRate = c.m_MagneticRate;
	const auto append = [&](CartesianSubgridConductivityData::Rows &rows, std::size_t edge, double value)
	{
		rows.edge.push_back(index(edge));
		rows.weight.push_back(scalar(value));
		rows.low.push_back(scalar(value - double(rows.weight.back())));
	};
	std::vector<std::vector<std::pair<std::size_t, double>>> rows(m_Grid.EdgeCount());
	for (std::size_t edge = 0; edge < rows.size(); ++edge)
	{
		if (c.m_Diagonal[edge] != 0)
		{
			rows[edge].emplace_back(edge, c.m_Diagonal[edge]);
		}
	}
	for (const auto &term : c.m_OffDiagonal)
	{
		rows[term.first.first].emplace_back(term.first.second, term.second);
		rows[term.first.second].emplace_back(term.first.first, term.second);
	}
	for (auto &row : rows)
	{
		std::sort(row.begin(), row.end());
		for (const auto &term : row)
		{
			append(loss.conductance, term.first, term.second);
		}
		loss.conductance.start.push_back(index(loss.conductance.edge.size()));
	}
	for (const auto &trace : c.m_Trace)
	{
		for (const auto &term : trace)
		{
			append(loss.trace, term.first, term.second);
		}
		loss.trace.start.push_back(index(loss.trace.edge.size()));
	}
	return GPU_CreateCartesianSubgridAdvanced(backend, m_Grid, std::move(data),
		GPU_CreateCartesianSubgridConductivity(std::move(loss)), c.m_Timestep);
}

void Engine_Ext_CartesianSubgrid_Conductivity::ReadDeviceDiagnostics(GPU_CartesianSubgrid &device)
{
	if (m_Conductivity->HasLoss())
	{
		const auto loss = device.ReadAdvanceLedger();
		m_ElectricLoss = {loss[0], 0};
		m_MagneticLoss = {loss[1], 0};
	}
}

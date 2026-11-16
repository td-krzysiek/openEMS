/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "engine_ext_cartesiansubgrid_conductivity.h"

#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>

Engine_Ext_CartesianSubgrid_PML::State::State(std::size_t size)
	: history(size), correction(size), previous(size), target(size)
{
}

Engine_Ext_CartesianSubgrid_PML::Engine_Ext_CartesianSubgrid_PML(Operator_Ext_CartesianSubgrid_PML *op)
	: Engine_Extension(op), m_PML(op), m_Electric(op->m_Electric.size()), m_Magnetic(op->m_Magnetic.size())
{
	SetPriority(ENG_EXT_PRIO_UPML);
	op->PreparedFootprint();
}

void Engine_Ext_CartesianSubgrid_PML::Prepare(
	const std::vector<Operator_Ext_CartesianSubgrid_PML::Derivative> &nodes, State &state,
	const std::vector<double> &source, const std::vector<double> &field, double sign)
{
	const double dt = m_PML->m_Compact->GetTimestep();
	for (std::size_t i = 0; i < nodes.size(); ++i)
	{
		const auto &node = nodes[i];
		double derivative = 0;
		for (const auto &term : node.terms)
			derivative += term.value * source[term.source];
		// The zero-frequency-shift, unit-kappa recursive convolution:
		// psi_new = exp(-rate*dt)*psi_old + (exp(-rate*dt)-1)*D(field).
		state.history[i] = node.decay * state.history[i] + node.change * derivative;
		state.correction[i] = sign * dt * state.history[i];
		state.previous[i] = field[node.target];
	}
}

void Engine_Ext_CartesianSubgrid_PML::Finish(
	const std::vector<Operator_Ext_CartesianSubgrid_PML::Derivative> &nodes, State &state,
	const std::vector<double> &field)
{
	// Several derivative directions may share a target. Each correction uses
	// the same natural field; the work ledger is evaluated in application order.
	for (std::size_t i = 0; i < nodes.size(); ++i)
	{
		state.target[i] = field[nodes[i].target];
		for (std::size_t j = i; j > 0 && nodes[j - 1].target == nodes[i].target; --j)
			state.target[i] += state.correction[j - 1];
		const double before = state.target[i];
		state.target[i] += state.correction[i];
		if (!std::isfinite(state.target[i]))
			throw std::runtime_error("CartesianSubgrid: nonfinite PML field");
		// After a prior correction, both endpoints in the electric/magnetic
		// work identity acquire that prior change.
		const double prior = before - field[nodes[i].target];
		m_WorkOut -= .5 * state.correction[i] * (state.previous[i] + prior + state.target[i]);
	}
}

void Engine_Ext_CartesianSubgrid_PML::Apply(
	const std::vector<Operator_Ext_CartesianSubgrid_PML::Derivative> &nodes, const State &state,
	std::vector<double> &field)
{
	for (std::size_t i = 0; i < nodes.size(); ++i)
		field[nodes[i].target] = state.target[i];
}

void Engine_Ext_CartesianSubgrid_PML::DoPreVoltageUpdates()
{
	m_HostAdvanced = true;
	const auto &g = m_PML->m_Compact->Grid();
	Prepare(m_PML->m_Electric, m_Electric, g.m_Magnetic, g.m_Electric, 1);
}
void Engine_Ext_CartesianSubgrid_PML::DoPostVoltageUpdates()
{
	Finish(m_PML->m_Electric, m_Electric, m_PML->m_Compact->Grid().m_Electric);
}
void Engine_Ext_CartesianSubgrid_PML::Apply2Voltages()
{
	Apply(m_PML->m_Electric, m_Electric, m_PML->m_Compact->Grid().m_Electric);
}
void Engine_Ext_CartesianSubgrid_PML::DoPreCurrentUpdates()
{
	m_HostAdvanced = true;
	const auto &g = m_PML->m_Compact->Grid();
	Prepare(m_PML->m_Magnetic, m_Magnetic, g.m_Electric, g.m_Magnetic, -1);
}
void Engine_Ext_CartesianSubgrid_PML::DoPostCurrentUpdates()
{
	Finish(m_PML->m_Magnetic, m_Magnetic, m_PML->m_Compact->Grid().m_Magnetic);
}
void Engine_Ext_CartesianSubgrid_PML::Apply2Current()
{
	Apply(m_PML->m_Magnetic, m_Magnetic, m_PML->m_Compact->Grid().m_Magnetic);
}

CartesianSubgridPMLData Engine_Ext_CartesianSubgrid_PML::PrepareDeviceData() const
{
	if (m_HostAdvanced || m_DeviceConfigured || m_WorkOut != 0 ||
		std::any_of(m_Electric.history.begin(), m_Electric.history.end(), [](double value) { return value != 0; }) ||
		std::any_of(m_Magnetic.history.begin(), m_Magnetic.history.end(), [](double value) { return value != 0; }))
	{
		throw std::invalid_argument("CartesianSubgrid: device PML creation requires fresh paired state");
	}
	const auto index = [](std::size_t value)
	{
		if (value > std::numeric_limits<std::uint32_t>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: device PML table exceeds 32-bit indices");
		}
		return static_cast<std::uint32_t>(value);
	};
	const auto scalar = [](double value)
	{
		const float converted = float(value);
		if (!std::isfinite(value) || !std::isfinite(converted) ||
			(value != 0 && std::abs(converted) < std::numeric_limits<float>::min()))
		{
			throw std::invalid_argument("CartesianSubgrid: PML coefficient is not normal FP32");
		}
		return converted;
	};
	CartesianSubgridPMLData data;
	data.timestep = m_PML->m_Compact->GetTimestep();
	const auto rows = [&](const auto &nodes, CartesianSubgridPMLData::Rows &output)
	{
		for (std::size_t row = 0; row < nodes.size(); ++row)
		{
			const auto &node = nodes[row];
			if (row && node.target != nodes[row - 1].target)
			{
				output.group.push_back(index(row));
			}
			output.target.push_back(index(node.target));
			output.decay.push_back(scalar(node.decay));
			output.change.push_back(scalar(node.change));
			for (const auto &term : node.terms)
			{
				const float high = scalar(term.value);
				output.source.push_back(index(term.source));
				output.weight.push_back(high);
				output.low.push_back(scalar(term.value - double(high)));
			}
			output.start.push_back(index(output.source.size()));
		}
		if (!nodes.empty())
		{
			output.group.push_back(index(nodes.size()));
		}
	};
	rows(m_PML->m_Electric, data.electric);
	rows(m_PML->m_Magnetic, data.magnetic);
	return data;
}

const CartesianSubgrid &Engine_Ext_CartesianSubgrid_PML::Grid() const
{
	return m_PML->m_Compact->Grid();
}

CartesianSubgridPMLData Engine_Ext_CartesianSubgrid_PML::PrepareDomainData(
	const Engine_Ext_CartesianSubgrid_Conductivity *advance) const
{
	if (advance && &advance->Grid() != &Grid())
	{
		throw std::invalid_argument("CartesianSubgrid: PML and domain advance have different grids");
	}
	for (std::size_t k = 0; k < m_PML->m_Compact->GetNumberOfExtentions(); ++k)
	{
		const auto *extension = m_PML->m_Compact->GetExtension(k);
		if (extension != m_PML && (!advance || extension != &advance->PreparedOperator()))
		{
			throw std::invalid_argument("CartesianSubgrid: other PML domain histories are not prepared");
		}
	}
	if (advance)
	{
		advance->PreparedOperator().ValidatePreparedExtensions();
	}
	return PrepareDeviceData();
}

std::unique_ptr<GPU_CartesianSubgridPML> Engine_Ext_CartesianSubgrid_PML::CreateDevice(
	GPU_Backend &backend, const GPU_CartesianSubgridFields &fields)
{
	auto device = GPU_CreateCartesianSubgridPML(backend, fields, PrepareDeviceData());
	m_DeviceConfigured = true;
	return device;
}

void Engine_Ext_CartesianSubgrid_PML::ReadDeviceDiagnostics(GPU_CartesianSubgridPML &device)
{
	if (!m_DeviceConfigured)
	{
		throw std::invalid_argument("CartesianSubgrid: no paired device PML state");
	}
	m_WorkOut = device.ReadBoundaryWorkOut();
}

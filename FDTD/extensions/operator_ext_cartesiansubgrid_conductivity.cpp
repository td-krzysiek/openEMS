/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "operator_ext_cartesiansubgrid_conductivity.h"
#include "engine_ext_cartesiansubgrid_conductivity.h"
#include "operator_ext_cartesiansubgrid_pml.h"
#include "operator_ext_cartesiansubgrid_mur.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/subgrid_interface_limits.h"
#include "FDTD/subgrid_setup_timing.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>
#include <limits>

namespace
{
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}
}

Operator_Ext_CartesianSubgrid_Conductivity::Operator_Ext_CartesianSubgrid_Conductivity(
	Operator *op, Sampler sampler)
	: Operator_Ext_CartesianSubgrid_Advance(op), m_Compact(dynamic_cast<Operator_CartesianSubgrid *>(op)),
	  m_Sample(std::move(sampler))
{
	Require(m_Compact && bool(m_Sample), "CartesianSubgrid conductivity requires a compact operator and sampler");
}

void Operator_Ext_CartesianSubgrid_Conductivity::SetSurfaceAdmittances(
	std::vector<SurfaceAdmittance> surfaces)
{
	Require(m_Timestep == 0, "CartesianSubgrid surface admittance must be configured before preparation");
	m_Surfaces = std::move(surfaces);
}

bool Operator_Ext_CartesianSubgrid_Conductivity::BuildExtension()
{
	CartesianSubgridSetupTiming timing("conductivity");
	const auto &g = m_Compact->Grid();
	m_Timestep = m_Compact->GetTimestep();
	g.CheckStep(m_Timestep);
	Require(g.m_PendingMagneticTimestep == 0, "CartesianSubgrid conductivity setup requires a complete timestep");
	Require(g.m_PortVoltage.empty() || g.m_PortTimestep == m_Timestep,
		"CartesianSubgrid conductivity/port timesteps differ");
	m_HasLoss = false;
	m_Diagonal.assign(g.EdgeCount(), 0);
	m_MagneticRate.assign(g.FaceCount(), 0);
	m_IndependentInverse.assign(g.EdgeCount(), 1);
	m_OffDiagonal.clear();
	m_Trace.clear();
	m_Factors.clear();
	m_Responses.clear();
	m_PortLower.clear();
	m_SurfaceTrace.clear();
	m_SurfaceRate.clear();
	m_SurfaceDenominator.clear();
	timing.Phase("coefficient_allocation");
	for (std::size_t c = 0; c < g.CellCount(); ++c)
	{
		if (g.LeafIsOpaque(c))
		{
			continue;
		}
		const auto start = g.LeafStart(c), length = g.LeafLengths(c);
		auto stop = start;
		for (unsigned int n = 0; n < 3; ++n)
		{
			stop[n] += length[n];
		}
		const auto sigma = m_Sample(start, stop);
		Require(sigma.electric >= 0 && std::isfinite(sigma.electric) &&
			sigma.magnetic >= 0 && std::isfinite(sigma.magnetic),
			"CartesianSubgrid requires finite nonnegative conductivities");
		Require((!sigma.electricDiagonal || sigma.Electric(0) == sigma.electric) &&
			(!sigma.magneticDiagonal || sigma.Magnetic(0) == sigma.magnetic),
			"CartesianSubgrid diagonal conductivity component zero must match its scalar alias");
		for (unsigned int n = 0; n < 3; ++n)
		{
			Require(sigma.Electric(n) >= 0 && std::isfinite(sigma.Electric(n)) &&
				sigma.Magnetic(n) >= 0 && std::isfinite(sigma.Magnetic(n)),
				"CartesianSubgrid requires finite nonnegative diagonal conductivities");
			m_HasLoss = m_HasLoss || sigma.Electric(n) != 0 || sigma.Magnetic(n) != 0;
			const double rate = sigma.Magnetic(n) / g.m_Media[c].Permeability(n);
			Require(std::isfinite(rate) && std::isfinite(.5 * m_Timestep * rate),
				"CartesianSubgrid magnetic conductivity rate overflows");
			m_MagneticRate[6 * c + 2 * n] = rate;
			m_MagneticRate[6 * c + 2 * n + 1] = rate;
		}
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (sigma.Electric(n) == 0)
			{
				continue;
			}
			const double weight = std::sqrt(sigma.Electric(n) * length[0] * length[1] * length[2] / 4);
			Require(weight > 0 && std::isfinite(weight), "CartesianSubgrid conductivity trace weight is unrepresentable");
			for (unsigned int i = 0; i < 2; ++i)
			{
				for (unsigned int j = 0; j < 2; ++j)
				{
					const auto local = g.Restrict(g.LocalEdge(g.m_Leaves[c], n, i, j));
					Trace trace;
					for (const auto &term : local.weights)
					{
						const double value = weight * term.second / g.m_SqrtMass[term.first];
						Require(std::isfinite(value) && value != 0, "CartesianSubgrid normalized conductivity trace is unrepresentable");
						trace.emplace_back(term.first, value);
					}
					for (const auto &a : trace)
					{
						m_Diagonal[a.first] += a.second * a.second;
						for (const auto &b : trace)
						{
							if (a.first < b.first)
							{
								m_OffDiagonal[{a.first, b.first}] += a.second * b.second;
							}
						}
					}
					if (!trace.empty())
					{
						m_Trace.push_back(std::move(trace));
					}
				}
			}
		}
	}
	timing.Phase("volume_traces");
	// Physical conductance and eliminated inductive states have distinct ledgers.
	const auto addTrace = [](std::vector<double> &diagonal, Couplings &off,
		const auto &trace, double weight)
	{
		for (const auto &a : trace)
		{
			diagonal[a.first] += weight * a.second * a.second;
			for (const auto &b : trace)
			{
				if (a.first < b.first)
				{
					off[{a.first, b.first}] += weight * a.second * b.second;
				}
			}
		}
	};
	for (const auto &surface : m_Surfaces)
	{
		Require(std::isfinite(surface.conductance) && surface.conductance >= 0,
			"CartesianSubgrid surface conductance must be finite and nonnegative");
		CartesianSubgrid::Integral path;
		if (surface.localEdge)
		{
			const auto edge = *surface.localEdge;
			path = g.LeafEdgeVoltage(edge.leaf, edge.axis, edge.first, edge.second);
			const auto origin = g.LeafStart(edge.leaf), length = g.LeafLengths(edge.leaf);
			auto start = origin;
			start[(edge.axis + 1) % 3] += edge.first * length[(edge.axis + 1) % 3];
			start[(edge.axis + 2) % 3] += edge.second * length[(edge.axis + 2) % 3];
			auto stop = start;
			stop[edge.axis] += length[edge.axis];
			for (unsigned int n = 0; n < 3; ++n)
			{
				const double tolerance = std::min(length[n] * 1e-8,
					32 * std::numeric_limits<double>::epsilon() *
						std::max({std::abs(start[n]), std::abs(stop[n]), length[n]}));
				Require(std::isfinite(surface.start[n]) && std::isfinite(surface.stop[n]) &&
					std::abs(surface.start[n] - start[n]) <= tolerance &&
					std::abs(surface.stop[n] - stop[n]) <= tolerance,
					"CartesianSubgrid surface edge index/physical coordinates differ");
			}
		}
		else
		{
			path = g.VoltagePath(surface.start, surface.stop);
		}
		Require(!path.weights.empty(), "CartesianSubgrid surface voltage path is empty or PEC");
		if (surface.conductance > 0)
		{
			Trace trace = path.weights;
			for (auto &term : trace)
			{
				term.second *= std::sqrt(surface.conductance);
			}
			addTrace(m_Diagonal, m_OffDiagonal, trace, 1);
			m_Trace.push_back(std::move(trace));
			m_HasLoss = true;
		}
		for (const auto &pole : surface.poles)
		{
			Require(std::isfinite(pole.resistance) && pole.resistance >= 0 &&
				std::isfinite(pole.inductance) && pole.inductance > 0,
				"CartesianSubgrid surface pole requires finite R>=0 and L>0");
			Trace trace = path.weights;
			for (auto &term : trace)
			{
				term.second /= std::sqrt(pole.inductance);
				Require(std::isfinite(term.second) && term.second != 0,
					"CartesianSubgrid surface trace is unrepresentable");
			}
			const double rate = pole.resistance / pole.inductance;
			const double denominator = 1 + .5 * m_Timestep * rate;
			Require(std::isfinite(rate) && std::isfinite(denominator),
				"CartesianSubgrid surface relaxation overflows");
			m_SurfaceTrace.push_back(std::move(trace));
			m_SurfaceRate.push_back(rate);
			m_SurfaceDenominator.push_back(denominator);
			m_HasLoss = true; // Includes lossless stored inductive energy.
		}
	}
	timing.Phase("surface_traces");
	m_EffectiveDiagonal = m_Diagonal;
	m_EffectiveOffDiagonal = m_OffDiagonal;
	for (std::size_t pole = 0; pole < m_SurfaceTrace.size(); ++pole)
	{
		addTrace(m_EffectiveDiagonal, m_EffectiveOffDiagonal, m_SurfaceTrace[pole],
			.5 * m_Timestep / m_SurfaceDenominator[pole]);
	}
	if (!m_HasLoss)
	{
		timing.Complete();
		return true;
	}
	timing.Phase("effective_conductance");
	// Mur remains outside the qualified joint recurrence. PML compatibility is
	// checked after all operators have prepared their actual protected supports.
	for (std::size_t i = 0; i < m_Compact->GetNumberOfExtentions(); ++i)
	{
		auto *extension = m_Compact->GetExtension(i);
		Require(!dynamic_cast<Operator_Ext_CartesianSubgrid_Mur *>(extension),
			"CartesianSubgrid conductivity with Mur histories is not yet qualified");
	}
	Couplings mass;
	for (const auto &term : g.m_MassOffDiagonal)
	{
		mass.emplace(term.first, term.second);
	}
	for (const auto &term : m_EffectiveOffDiagonal)
	{
		if (!mass.count(term.first))
		{
			throw std::invalid_argument("CartesianSubgrid conductivity escaped the electric mass graph: edges " +
				std::to_string(term.first.first) + "," + std::to_string(term.first.second) +
				" weight " + std::to_string(term.second));
		}
	}
	for (std::size_t i = 0; i < g.EdgeCount(); ++i)
	{
		const double diagonal = 1 + .5 * m_Timestep * m_EffectiveDiagonal[i];
		Require(diagonal > 0 && std::isfinite(diagonal), "CartesianSubgrid implicit conductivity diagonal overflows");
		m_IndependentInverse[i] = 1 / diagonal;
	}
	for (const auto &original : g.m_Blocks)
	{
		auto block = original;
		const auto n = block.edges.size(), stride = block.halfBandwidth + 1;
		block.lower.assign(CartesianSubgridDetail::InterfaceFactorEntries(n, block.halfBandwidth), 0);
		for (std::size_t i = 0; i < n; ++i)
		{
			const auto edge = block.edges[i];
			m_IndependentInverse[edge] = 1;
			block.lower[i * stride] = 1 + .5 * m_Timestep * m_EffectiveDiagonal[edge];
			for (std::size_t j = i > block.halfBandwidth ? i - block.halfBandwidth : 0; j < i; ++j)
			{
				const auto pair = std::minmax(edge, block.edges[j]);
				const std::pair<std::size_t, std::size_t> key{pair.first, pair.second};
				const auto a = mass.find(key), b = m_EffectiveOffDiagonal.find(key);
				block.lower[i * stride + i - j] = (a != mass.end() ? a->second : 0) +
					.5 * m_Timestep * (b != m_EffectiveOffDiagonal.end() ? b->second : 0);
			}
		}
		Factor(block);
		m_Factors.push_back(std::move(block));
	}
	timing.Phase("implicit_factors");
	BuildPortResponses();
	timing.Phase("port_responses");
	timing.Complete();
	return true;
}

void Operator_Ext_CartesianSubgrid_Conductivity::ValidatePreparedExtensions() const
{
	if (!m_HasLoss)
	{
		return;
	}
	const auto &grid = m_Compact->Grid();
	std::vector<bool> affected(grid.EdgeCount(), false);
	for (std::size_t edge = 0; edge < affected.size(); ++edge)
	{
		affected[edge] = m_EffectiveDiagonal[edge] != 0;
	}
	// Any nonzero G on a coupled mass block protects the entire block: an
	// implicit solve can change every member even if its own G diagonal is zero.
	for (const auto &block : grid.m_Blocks)
	{
		const bool coupled = std::any_of(block.edges.begin(), block.edges.end(),
			[&](std::size_t edge) { return affected[edge]; });
		if (coupled)
		{
			for (std::size_t edge : block.edges)
			{
				affected[edge] = true;
			}
		}
	}
	for (std::size_t i = 0; i < m_Compact->GetNumberOfExtentions(); ++i)
	{
		const auto *pml = dynamic_cast<Operator_Ext_CartesianSubgrid_PML *>(m_Compact->GetExtension(i));
		if (!pml)
		{
			continue;
		}
		const auto &footprint = pml->PreparedFootprint();
		for (std::size_t edge : footprint.electric)
		{
			Require(edge < affected.size() && !affected[edge],
				"CartesianSubgrid: conductivity/PML electric mass support overlap");
		}
		for (std::size_t face : footprint.magnetic)
		{
			Require(face < m_MagneticRate.size() && m_MagneticRate[face] == 0,
				"CartesianSubgrid: conductivity/PML magnetic support overlap");
		}
	}
}

void Operator_Ext_CartesianSubgrid_Conductivity::Factor(CartesianSubgrid::MassBlock &block)
{
	const auto n = block.edges.size(), band = block.halfBandwidth, stride = band + 1;
	for (std::size_t i = 0; i < n; ++i)
	{
		const auto first = i > band ? i - band : 0;
		for (std::size_t j = first; j <= i; ++j)
		{
			double value = block.lower[i * stride + i - j];
			for (std::size_t k = first; k < j; ++k)
			{
				value -= block.lower[i * stride + i - k] * block.lower[j * stride + j - k];
			}
			Require(std::isfinite(value) && (i != j || value > 0), "CartesianSubgrid conductivity factor is not SPD");
			block.lower[i * stride + i - j] = i == j ? std::sqrt(value) : value / block.lower[j * stride];
		}
	}
}

void Operator_Ext_CartesianSubgrid_Conductivity::Solve(std::vector<double> &values, std::vector<double> &scratch) const
{
	for (std::size_t i = 0; i < values.size(); ++i)
	{
		values[i] *= m_IndependentInverse[i];
	}
	for (const auto &block : m_Factors)
	{
		scratch.resize(block.edges.size());
		for (std::size_t i = 0; i < block.edges.size(); ++i)
		{
			scratch[i] = values[block.edges[i]];
		}
		block.Solve(scratch);
		for (std::size_t i = 0; i < block.edges.size(); ++i)
		{
			values[block.edges[i]] = scratch[i];
		}
	}
}

void Operator_Ext_CartesianSubgrid_Conductivity::BuildPortResponses()
{
	const auto &g = m_Compact->Grid();
	std::vector<double> scratch;
	for (std::size_t p = 0; p < g.m_PortVoltage.size(); ++p)
	{
		std::vector<double> response(g.EdgeCount(), 0);
		for (const auto &term : g.m_PortVoltage[p].weights)
		{
			response[term.first] = std::sqrt(g.m_PortConductance[p]) * term.second;
		}
		Solve(response, scratch);
		CartesianSubgrid::Integral sparse;
		for (std::size_t i = 0; i < response.size(); ++i)
		{
			if (response[i] != 0)
			{
				sparse.weights.emplace_back(i, response[i]);
			}
		}
		m_Responses.push_back(std::move(sparse));
	}
	const auto count = m_Responses.size();
	m_PortLower.assign(count * count, 0);
	for (std::size_t i = 0; i < count; ++i)
	{
		for (std::size_t j = 0; j <= i; ++j)
		{
			double value = i == j ? 1 : 0;
			for (const auto &a : g.m_PortVoltage[i].weights)
			{
				for (const auto &b : m_Responses[j].weights)
				{
					if (a.first == b.first)
					{
						value += .5 * m_Timestep * std::sqrt(g.m_PortConductance[i]) * a.second * b.second;
					}
				}
			}
			for (std::size_t k = 0; k < j; ++k)
			{
				value -= m_PortLower[i * count + k] * m_PortLower[j * count + k];
			}
			Require(std::isfinite(value) && (i != j || value > 0), "CartesianSubgrid conductivity port factor is not SPD");
			m_PortLower[i * count + j] = i == j ? std::sqrt(value) : value / m_PortLower[j * count + j];
		}
	}
}

Engine_Ext_CartesianSubgrid_Advance *Operator_Ext_CartesianSubgrid_Conductivity::CreateEngineExtention()
{
	m_Eng_Ext = new Engine_Ext_CartesianSubgrid_Conductivity(this);
	return static_cast<Engine_Ext_CartesianSubgrid_Advance *>(m_Eng_Ext);
}

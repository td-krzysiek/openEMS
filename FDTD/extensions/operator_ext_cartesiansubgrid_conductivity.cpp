/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "operator_ext_cartesiansubgrid_conductivity.h"
#include "engine_ext_cartesiansubgrid_conductivity.h"
#include "operator_ext_cartesiansubgrid_pml.h"
#include "operator_ext_cartesiansubgrid_mur.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/subgrid_interface_limits.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

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

bool Operator_Ext_CartesianSubgrid_Conductivity::BuildExtension()
{
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
	for (std::size_t c = 0; c < g.CellCount(); ++c)
	{
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
		m_HasLoss = m_HasLoss || sigma.electric != 0 || sigma.magnetic != 0;
		const double rate = sigma.magnetic / g.m_Media[c].permeability;
		Require(std::isfinite(rate) && std::isfinite(.5 * m_Timestep * rate), "CartesianSubgrid magnetic conductivity rate overflows");
		std::fill(m_MagneticRate.begin() + 6 * c, m_MagneticRate.begin() + 6 * (c + 1), rate);
		if (sigma.electric == 0)
		{
			continue;
		}
		const double weight = std::sqrt(sigma.electric * length[0] * length[1] * length[2] / 4);
		Require(weight > 0 && std::isfinite(weight), "CartesianSubgrid conductivity trace weight is unrepresentable");
		for (unsigned int n = 0; n < 3; ++n)
		{
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
	if (!m_HasLoss)
	{
		return true;
	}
	// Histories need a joint derivation before these combinations are admitted.
	for (std::size_t i = 0; i < m_Compact->GetNumberOfExtentions(); ++i)
	{
		auto *extension = m_Compact->GetExtension(i);
		Require(!dynamic_cast<Operator_Ext_CartesianSubgrid_PML *>(extension) &&
			!dynamic_cast<Operator_Ext_CartesianSubgrid_Mur *>(extension),
			"CartesianSubgrid conductivity with absorber histories is not yet qualified");
	}
	Couplings mass;
	for (const auto &term : g.m_MassOffDiagonal)
	{
		mass.emplace(term.first, term.second);
	}
	for (const auto &term : m_OffDiagonal)
	{
		Require(mass.count(term.first), "CartesianSubgrid conductivity escaped the electric mass graph");
	}
	for (std::size_t i = 0; i < g.EdgeCount(); ++i)
	{
		const double diagonal = 1 + .5 * m_Timestep * m_Diagonal[i];
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
			block.lower[i * stride] = 1 + .5 * m_Timestep * m_Diagonal[edge];
			for (std::size_t j = i > block.halfBandwidth ? i - block.halfBandwidth : 0; j < i; ++j)
			{
				const auto pair = std::minmax(edge, block.edges[j]);
				const std::pair<std::size_t, std::size_t> key{pair.first, pair.second};
				const auto a = mass.find(key), b = m_OffDiagonal.find(key);
				block.lower[i * stride + i - j] = (a != mass.end() ? a->second : 0) +
					.5 * m_Timestep * (b != m_OffDiagonal.end() ? b->second : 0);
			}
		}
		Factor(block);
		m_Factors.push_back(std::move(block));
	}
	BuildPortResponses();
	return true;
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

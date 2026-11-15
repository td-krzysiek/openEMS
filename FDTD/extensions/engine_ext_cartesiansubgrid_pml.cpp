/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/operator_cartesiansubgrid.h"

#include <cmath>
#include <stdexcept>

Engine_Ext_CartesianSubgrid_PML::State::State(std::size_t size)
	: history(size), correction(size), previous(size), target(size)
{
}

Engine_Ext_CartesianSubgrid_PML::Engine_Ext_CartesianSubgrid_PML(Operator_Ext_CartesianSubgrid_PML *op)
	: Engine_Extension(op), m_PML(op), m_Electric(op->m_Electric.size()), m_Magnetic(op->m_Magnetic.size())
{
	SetPriority(ENG_EXT_PRIO_UPML);
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

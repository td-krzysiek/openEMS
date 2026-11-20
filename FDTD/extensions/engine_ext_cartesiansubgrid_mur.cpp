/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "engine_ext_cartesiansubgrid_mur.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "operator_ext_cartesiansubgrid_mur.h"

#include <cmath>
#include <stdexcept>

Engine_Ext_CartesianSubgrid_Mur::Engine_Ext_CartesianSubgrid_Mur(Operator_Ext_CartesianSubgrid_Mur *op)
	: Engine_Extension(op), m_Mur(op), m_Previous(op->m_Nodes.size()), m_Target(op->m_Nodes.size()),
	  m_Natural(op->m_Nodes.size())
{
}

void Engine_Ext_CartesianSubgrid_Mur::DoPreVoltageUpdates()
{
	const auto &grid = m_Mur->m_Compact->Grid();
	for (std::size_t i = 0; i < m_Target.size(); ++i)
	{
		const auto &node = m_Mur->m_Nodes[i];
		m_Previous[i] = grid.Evaluate(node.boundary);
		m_Target[i] = grid.Evaluate(node.inside) - m_Mur->m_Coefficient * m_Previous[i];
	}
}

void Engine_Ext_CartesianSubgrid_Mur::DoPostVoltageUpdates()
{
	const auto &grid = m_Mur->m_Compact->Grid();
	for (std::size_t i = 0; i < m_Target.size(); ++i)
	{
		m_Target[i] += m_Mur->m_Coefficient * grid.Evaluate(m_Mur->m_Nodes[i].inside);
		m_Natural[i] = grid.Evaluate(m_Mur->m_Nodes[i].boundary);
	}
}

void Engine_Ext_CartesianSubgrid_Mur::Apply2Voltages()
{
	auto &grid = m_Mur->m_Compact->Grid();
	for (std::size_t i = 0; i < m_Target.size(); ++i)
	{
		if (!std::isfinite(m_Target[i]))
		{
			throw std::runtime_error("CartesianSubgrid: nonfinite Mur field");
		}
		const auto &node = m_Mur->m_Nodes[i];
		const double delta = m_Target[i] - m_Natural[i];
		// The two-cell collar has diagonal mass and no source/circuit work.
		// At the E/H seam, the leapfrog invariant's electric-correction work
		// is 1/2 * (V_corrected-V_natural)*(V_old+V_corrected)/weight^2.
		// This is signed boundary work; Mur is not per-step dissipative.
		m_WorkOut -= .5 * delta * (m_Previous[i] + m_Target[i]) * node.inverseWeightSquared;
		grid.AddElectric(node.edge, delta / node.length);
	}
}

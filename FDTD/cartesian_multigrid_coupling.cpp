/*
 * Derived from the cylindrical coupling in engine_subgrid.cpp.
 * Copyright (C) 2010-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
 * Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "cartesian_multigrid_coupling.h"
#include "engine.h"

#include <algorithm>
#include <stdexcept>

CartesianMultiGridCoupling::CartesianMultiGridCoupling(Engine *base, Engine *coarse,
													   const MultiGridInterpolation &interpolation,
													   const unsigned int lines[3])
	: m_Base(base), m_Coarse(coarse), m_Interpolation(interpolation)
{
	if (!base || !coarse || base->GetType() != coarse->GetType() || base->GetType() == Engine::GPU)
	{
		throw std::invalid_argument("Cartesian multi-grid coupling requires matching host engines");
	}
	std::copy(lines, lines + 3, m_Lines);
}

void CartesianMultiGridCoupling::SyncVoltages()
{
	const unsigned int x = m_Interpolation.split_pos - 1;
	for (unsigned int y = 0; y + 1 < m_Lines[1]; y += 2)
	{
		for (unsigned int z = 0; z < m_Lines[2]; ++z)
		{
			m_Coarse->SetVolt(0, x, y / 2, z, 0);
			m_Coarse->SetVolt(2, x, y / 2, z, m_Base->GetVolt(2, x, y, z));
			m_Coarse->SetVolt(1, x, y / 2, z, m_Base->GetVolt(1, x, y, z) + m_Base->GetVolt(1, x, y + 1, z));
		}
	}
}

void CartesianMultiGridCoupling::SyncCurrents()
{
	// R maps fine boundary voltages to coarse: Ez is sampled at even Y,
	// Ey is the sum of its two fine edges. Transfer the conjugate currents
	// by R^T, so Vf^T If == (R Vf)^T Ic. Point interpolation of these
	// interface currents is not the adjoint and can inject energy.
	const unsigned int x = m_Interpolation.split_pos - 2;
	for (unsigned int y = 0; y < m_Lines[1]; ++y)
	{
		for (unsigned int z = 0; z < m_Lines[2]; ++z)
		{
			m_Base->SetCurr(1, x, y, z, y % 2 == 0 ? m_Coarse->GetCurr(1, x, y / 2, z) : 0);
			m_Base->SetCurr(2, x, y, z, y + 1 < m_Lines[1] ? m_Coarse->GetCurr(2, x, y / 2, z) : 0);
		}
	}
}

void CartesianMultiGridCoupling::Interpolate()
{
	for (unsigned int x = 0; x + 1 < m_Interpolation.split_pos; ++x)
	{
		InterpolatePlane(x, false);
	}
	for (unsigned int x = 0; x + 2 < m_Interpolation.split_pos; ++x)
	{
		InterpolatePlane(x, true);
	}
}

void CartesianMultiGridCoupling::InterpolatePlane(unsigned int x, bool currents)
{
	const auto &ip = m_Interpolation;
	const auto *first = currents ? ip.pos_i_2p : ip.pos_v_2p;
	const auto *second = currents ? ip.pos_i_2pp : ip.pos_v_2pp;
	const auto *wFirst = currents ? ip.w_i_2p : ip.w_v_2p;
	const auto *wSecond = currents ? ip.w_i_2pp : ip.w_v_2pp;
	for (unsigned int y = 0; y < m_Lines[1]; ++y)
	{
		for (unsigned int z = 0; z < m_Lines[2]; ++z)
		{
			for (unsigned int component = 0; component < 3; ++component)
			{
				const unsigned int direction = component == 1 ? 1 : 0;
				const unsigned int a = first[direction][y];
				const unsigned int b = second[direction][y];
				const FDTD_FLOAT value = currents ? wFirst[direction][y] * m_Coarse->GetCurr(component, x, a, z) +
														wSecond[direction][y] * m_Coarse->GetCurr(component, x, b, z)
												  : wFirst[direction][y] * m_Coarse->GetVolt(component, x, a, z) +
														wSecond[direction][y] * m_Coarse->GetVolt(component, x, b, z);
				if (currents)
				{
					m_Base->SetCurr(component, x, y, z, value);
				}
				else
				{
					m_Base->SetVolt(component, x, y, z, value);
				}
			}
		}
	}
}

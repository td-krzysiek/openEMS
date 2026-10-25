/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "engine_ext_dispersive.h"
#include "operator_ext_dispersive.h"
#include "FDTD/engine_sse.h"

Engine_Ext_Dispersive::Engine_Ext_Dispersive(Operator_Ext_Dispersive* op_ext_disp) : Engine_Extension(op_ext_disp)
{
	m_Op_Ext_Disp = op_ext_disp;
	m_TilingSupported = true;
	m_Order = m_Op_Ext_Disp->m_Order;
	int order = m_Order;
	curr_ADE = new FDTD_FLOAT**[order];
	volt_ADE = new FDTD_FLOAT**[order];
	for (int o=0;o<order;++o)
	{
		curr_ADE[o] = new FDTD_FLOAT*[3];
		volt_ADE[o] = new FDTD_FLOAT*[3];
		for (int n=0; n<3; ++n)
		{
			if (m_Op_Ext_Disp->m_curr_ADE_On[o]==true)
			{
				curr_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count[o]];
				for (size_t i=0; i<m_Op_Ext_Disp->m_LM_Count[o]; ++i)
					curr_ADE[o][n][i]=0.0;
			}
			else
				curr_ADE[o][n] = NULL;
			if (m_Op_Ext_Disp->m_volt_ADE_On[o]==true)
			{
				volt_ADE[o][n] = new FDTD_FLOAT[m_Op_Ext_Disp->m_LM_Count[o]];
				for (size_t i=0; i<m_Op_Ext_Disp->m_LM_Count[o]; ++i)
					volt_ADE[o][n][i]=0.0;
			}
			else
				volt_ADE[o][n] = NULL;
		}
	}
}

Engine_Ext_Dispersive::~Engine_Ext_Dispersive()
{
	if (curr_ADE==NULL && volt_ADE==NULL)
		return;

	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		for (int n=0; n<3; ++n)
		{
			delete[] curr_ADE[o][n];
			delete[] volt_ADE[o][n];
		}
		delete[] curr_ADE[o];
		delete[] volt_ADE[o];
	}
	delete[] curr_ADE;
	curr_ADE=NULL;

	delete[] volt_ADE;
	volt_ADE=NULL;
}

template <typename EngType>
void Engine_Ext_Dispersive::Apply2VoltagesImpl(EngType* eng)
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_volt_ADE_On[o]==false) continue;

		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		for (size_t i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
		{
			eng->EngType::SetVolt(0,pos[0][i],pos[1][i],pos[2][i],
				eng->EngType::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][0][i]
			);
			eng->EngType::SetVolt(1,pos[0][i],pos[1][i],pos[2][i],
				eng->EngType::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][1][i]
			);
			eng->EngType::SetVolt(2,pos[0][i],pos[1][i],pos[2][i],
				eng->EngType::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][2][i]
			);
		}
	}
}

void Engine_Ext_Dispersive::Apply2Voltages()
{
	ENG_DISPATCH(Apply2VoltagesImpl);
}

template <typename EngType>
void Engine_Ext_Dispersive::Apply2CurrentImpl(EngType* eng)
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_curr_ADE_On[o]==false) continue;

		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		for (size_t i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
		{
			eng->EngType::SetCurr(0,pos[0][i],pos[1][i],pos[2][i],
				eng->EngType::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][0][i]
			);
			eng->EngType::SetCurr(1,pos[0][i],pos[1][i],pos[2][i],
				eng->EngType::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][1][i]
			);
			eng->EngType::SetCurr(2,pos[0][i],pos[1][i],pos[2][i],
				eng->EngType::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][2][i]
			);
		}
	}
}

void Engine_Ext_Dispersive::Apply2Current()
{
	ENG_DISPATCH(Apply2CurrentImpl);
}

void Engine_Ext_Dispersive::InitializeTiling(const std::vector<Range3D>& tiles)
{
	m_volt_map.clear();
	m_curr_map.clear();


	for (auto& tile : tiles)
	{
		for (int o = 0; o < m_Order; ++o)
		{
			if (!tile.HasVoltage() || !m_Op_Ext_Disp->m_volt_ADE_On[o])
				continue;

			unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];
			unsigned int start[3], stop[3];
			for (int n=0; n<3; ++n) { start[n]=tile.voltageStart[n]; stop[n]=tile.voltageStop[n]; }

			auto entry = m_volt_map.emplace(GetTileKey(o, start, stop), std::vector<int>());
			if (!entry.second) continue;
			auto& indices = entry.first->second;

			for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
			{
				if (InsideTile(start, stop, pos[0][i], pos[1][i], pos[2][i]))
				{
					indices.push_back(i);
				}
			}
		}

		for (int o = 0; o < m_Order; ++o)
		{
			if (!tile.HasCurrent() || !m_Op_Ext_Disp->m_curr_ADE_On[o])
				continue;

			unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];
			unsigned int start[3], stop[3];
			for (int n=0; n<3; ++n) { start[n]=tile.currentStart[n]; stop[n]=tile.currentStop[n]; }

			auto entry = m_curr_map.emplace(GetTileKey(o, start, stop), std::vector<int>());
			if (!entry.second) continue;
			auto& indices = entry.first->second;

			for (unsigned int i=0; i<m_Op_Ext_Disp->m_LM_Count.at(o); ++i)
			{
				if (InsideTile(start, stop, pos[0][i], pos[1][i], pos[2][i]))
				{
					indices.push_back(i);
				}
			}
		}
	}

}

bool Engine_Ext_Dispersive::InsideTile(
	const unsigned int start[3], const unsigned int stop[3],
	unsigned int ade_x, unsigned int ade_y, unsigned int ade_z
)
{
	if (ade_x < start[0] || ade_x > stop[0])
		return false;
	else if (ade_y < start[1] || ade_y > stop[1])
		return false;
	else if (ade_z < start[2] || ade_z > stop[2])
		return false;
	else
		return true;
}

void Engine_Ext_Dispersive::Apply2Voltages(int timestep, unsigned int start[3], unsigned int stop[3])
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_volt_ADE_On[o]==false) continue;

		const auto& vec = m_volt_map.at(GetTileKey(o, start, stop));
		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		Engine_sse* eng_sse = (Engine_sse*)m_Eng;
		for (auto& i : vec)
		{
			eng_sse->Engine_sse::SetVolt(0,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(0,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][0][i]);
			eng_sse->Engine_sse::SetVolt(1,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(1,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][1][i]);
			eng_sse->Engine_sse::SetVolt(2,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetVolt(2,pos[0][i],pos[1][i],pos[2][i]) - volt_ADE[o][2][i]);
		}
	}
}

void Engine_Ext_Dispersive::Apply2Current(int timestep, unsigned int start[3], unsigned int stop[3])
{
	for (int o=0;o<m_Op_Ext_Disp->m_Order;++o)
	{
		if (m_Op_Ext_Disp->m_curr_ADE_On[o]==false) continue;

		const auto& vec = m_curr_map.at(GetTileKey(o, start, stop));
		unsigned int **pos = m_Op_Ext_Disp->m_LM_pos[o];

		Engine_sse* eng_sse = (Engine_sse*)m_Eng;
		for (auto& i : vec)
		{
			eng_sse->Engine_sse::SetCurr(0,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(0,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][0][i]);
			eng_sse->Engine_sse::SetCurr(1,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(1,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][1][i]);
			eng_sse->Engine_sse::SetCurr(2,pos[0][i],pos[1][i],pos[2][i], eng_sse->Engine_sse::GetCurr(2,pos[0][i],pos[1][i],pos[2][i]) - curr_ADE[o][2][i]);
		}
	}
}

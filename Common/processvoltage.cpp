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

#include "processvoltage.h"
#include "FDTD/engine_interface_fdtd.h"
#include <iomanip>
#include <sstream>
#include "FDTD/operator_cartesiansubgrid.h"

ProcessVoltage::ProcessVoltage(Engine_Interface_Base* eng_if) : ProcessIntegral(eng_if)
{
}

ProcessVoltage::~ProcessVoltage()
{
}

std::string ProcessVoltage::GetIntegralName(int row) const
{
	if (row==0)
		return "voltage";
	return "unknown";
}

void ProcessVoltage::DefineStartStopCoord(double* dstart, double* dstop)
{
	m_LocalIntegral = {};
	if (auto *local = dynamic_cast<const Operator_CartesianSubgrid *>(Op))
	{
		const auto path = local->Grid().VoltagePath(local->DrawingPoint(dstart), local->DrawingPoint(dstop));
		m_LocalIntegral = [local, path]()
		{
			return local->Grid().Evaluate(path);
		};
		std::ostringstream description;
		description << std::setprecision(17)
					<< "% subgrid physical integral (m); base-index headers above are informational: ";
		const auto lo = local->DrawingPoint(dstart), hi = local->DrawingPoint(dstop);
		for (unsigned int n = 0; n < 3; ++n)
		{
			description << (n ? "," : "(") << lo[n] + local->GetDiscLine(n, 0) * local->GetGridDelta();
		}
		description << ") -> (";
		for (unsigned int n = 0; n < 3; ++n)
		{
			description << (n ? "," : "") << hi[n] + local->GetDiscLine(n, 0) * local->GetGridDelta();
		}
		description << ")";
		m_LocalGeometry = description.str();
	}
	Op->SnapToMesh(dstart, start, m_dualMesh, false, m_start_inside);
	Op->SnapToMesh(dstop, stop, m_dualMesh, false, m_stop_inside);
}

double ProcessVoltage::CalcIntegral()
{
	if (m_LocalIntegral)
	{
		return m_LocalIntegral();
	}
	//integrate voltages from start to stop on a line
	return m_Eng_Interface->CalcVoltageIntegral(start,stop);
}

std::string ProcessVoltage::IntegrationGeometry() const
{
	return m_LocalIntegral ? m_LocalGeometry + "\n" : "";
}

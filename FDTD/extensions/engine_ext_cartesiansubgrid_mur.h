/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef ENGINE_EXT_CARTESIANSUBGRID_MUR_H
#define ENGINE_EXT_CARTESIANSUBGRID_MUR_H

#include "engine_extension.h"
#include <vector>

class Operator_Ext_CartesianSubgrid_Mur;

class Engine_Ext_CartesianSubgrid_Mur : public Engine_Extension
{
  public:
	explicit Engine_Ext_CartesianSubgrid_Mur(Operator_Ext_CartesianSubgrid_Mur *op);
	void DoPreVoltageUpdates() override;
	void DoPostVoltageUpdates() override;
	void Apply2Voltages() override;
	//! Signed work removed by the discrete boundary correction, not radiated power.
	double BoundaryWorkOut() const
	{
		return m_WorkOut;
	}

  private:
	Operator_Ext_CartesianSubgrid_Mur *m_Mur;
	std::vector<double> m_Previous, m_Target, m_Natural;
	double m_WorkOut = 0;
};
#endif

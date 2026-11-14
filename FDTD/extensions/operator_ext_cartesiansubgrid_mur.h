/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPERATOR_EXT_CARTESIANSUBGRID_MUR_H
#define OPERATOR_EXT_CARTESIANSUBGRID_MUR_H

#include "FDTD/cartesian_subgrid.h"
#include "operator_extension.h"

class Operator_CartesianSubgrid;
class Engine_Ext_CartesianSubgrid_Mur;

//! First-order Mur on an unrefined, homogeneous, source-free boundary collar.
class Operator_Ext_CartesianSubgrid_Mur : public Operator_Extension
{
	friend class Engine_Ext_CartesianSubgrid_Mur;

  public:
	Operator_Ext_CartesianSubgrid_Mur(Operator *op, unsigned int face);
	void SetPhaseVelocity(double velocity);
	bool BuildExtension() override;
	Engine_Extension *CreateEngineExtention() override;
	std::string GetExtensionName() const override
	{
		return "Cartesian subgrid Mur ABC";
	}
	std::size_t NodeCount() const
	{
		return m_Nodes.size();
	}

  private:
	struct Node
	{
		CartesianSubgrid::Integral boundary, inside;
		std::size_t edge;
		double length, inverseWeightSquared;
	};
	Operator_CartesianSubgrid *m_Compact;
	unsigned int m_Face;
	double m_Velocity = 0, m_Coefficient = 0;
	std::vector<Node> m_Nodes;
};
#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPERATOR_EXT_CARTESIANSUBGRID_PML_H
#define OPERATOR_EXT_CARTESIANSUBGRID_PML_H

#include "FDTD/cartesian_subgrid.h"
#include "operator_extension.h"

class Operator_CartesianSubgrid;
class Engine_Ext_CartesianSubgrid_PML;

//! Classical stretched-coordinate PML (kappa=1, alpha=0), serial reference only.
//! A homogeneous, unrefined, source/PEC-free collar includes one inner guard cell.
class Operator_Ext_CartesianSubgrid_PML : public Operator_Extension
{
	friend class Engine_Ext_CartesianSubgrid_PML;

  public:
	explicit Operator_Ext_CartesianSubgrid_PML(Operator *op, const std::array<unsigned int, 6> &cells);
	bool BuildExtension() override;
	Engine_Extension *CreateEngineExtention() override;
	std::string GetExtensionName() const override { return "Cartesian subgrid convolution PML"; }
	std::size_t ElectricHistoryCount() const { return m_Electric.size(); }
	std::size_t MagneticHistoryCount() const { return m_Magnetic.size(); }

  private:
	struct Term
	{
		std::size_t source;
		double value;
	};
	struct Derivative
	{
		std::size_t target;
		double decay, change;
		std::vector<Term> terms;
	};
	Operator_CartesianSubgrid *m_Compact;
	std::array<unsigned int, 6> m_Cells;
	std::vector<Derivative> m_Electric, m_Magnetic;
};
#endif

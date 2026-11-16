/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPERATOR_EXT_CARTESIANSUBGRID_CONDUCTIVITY_H
#define OPERATOR_EXT_CARTESIANSUBGRID_CONDUCTIVITY_H

#include "cartesiansubgrid_advance.h"
#include "FDTD/cartesian_subgrid.h"
#include <functional>
#include <map>
#include <vector>

class Operator_CartesianSubgrid;
class Engine_Ext_CartesianSubgrid_Conductivity;

//! Constant, isotropic SI conductivities; not a broadband loss-tangent model.
class Operator_Ext_CartesianSubgrid_Conductivity : public Operator_Ext_CartesianSubgrid_Advance
{
	friend class Engine_Ext_CartesianSubgrid_Conductivity;
  public:
	struct Coefficients
	{
		double electric = 0; //!< S/m
		double magnetic = 0; //!< ohm/m, the magnetic conductivity of Maxwell's equations
	};
	using Sampler = std::function<Coefficients(const CartesianSubgrid::Point &, const CartesianSubgrid::Point &)>;
	Operator_Ext_CartesianSubgrid_Conductivity(Operator *op, Sampler sampler);
	bool BuildExtension() override;
	void ValidatePreparedExtensions() const override;
	Engine_Ext_CartesianSubgrid_Advance *CreateEngineExtention() override;
	std::string GetExtensionName() const override { return "Cartesian subgrid midpoint conductivity"; }
	bool HasLoss() const { return m_HasLoss; }

  private:
	using Couplings = std::map<std::pair<std::size_t, std::size_t>, double>;
	using Trace = std::vector<std::pair<std::size_t, double>>;
	void Solve(std::vector<double> &values, std::vector<double> &scratch) const;
	void BuildPortResponses();
	static void Factor(CartesianSubgrid::MassBlock &block);
	Operator_CartesianSubgrid *m_Compact;
	Sampler m_Sample;
	bool m_HasLoss = false;
	double m_Timestep = 0;
	std::vector<double> m_Diagonal, m_MagneticRate, m_IndependentInverse;
	Couplings m_OffDiagonal;
	std::vector<Trace> m_Trace;
	std::vector<CartesianSubgrid::MassBlock> m_Factors;
	std::vector<CartesianSubgrid::Integral> m_Responses;
	std::vector<double> m_PortLower;
};
#endif

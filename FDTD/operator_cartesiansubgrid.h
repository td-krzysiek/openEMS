/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPERATOR_CARTESIANSUBGRID_H
#define OPERATOR_CARTESIANSUBGRID_H

#include "cartesian_subgrid.h"
#include "operator.h"
#include <memory>

//! Native API adapter for the experimental compact XYZ numerical core.
class Operator_CartesianSubgrid : public Operator
{
  public:
	static Operator_CartesianSubgrid *New(const std::vector<double> &box, CartesianSubgrid::Index ratio = {2, 2, 2});
	void SetEngineType(Engine_Types::EngineType type) override;
	bool SetGeometryCSX(ContinuousStructure *geometry) override;
	void SetBoundaryCondition(int *conditions) override;
	void AddExtension(Operator_Extension *extension) override;
	int CalcECOperator(DebugFlags flags = None) override;
	Engine *CreateEngine() override;
	double GetNumberCells() const override;
	void ShowStat() const override;
	CartesianSubgrid &Grid() const
	{
		return *m_Grid;
	}
	CartesianSubgrid::Point Position(unsigned int n, const unsigned int pos[3], bool dual) const;
	Excitation *Signal() const
	{
		return m_Exc;
	}

  private:
	Operator_CartesianSubgrid(const std::vector<double> &box, CartesianSubgrid::Index ratio);
	void ValidateProperties() const;
	void PrepareSources();
	std::vector<double> m_Box;
	CartesianSubgrid::Index m_Ratio;
	std::unique_ptr<CartesianSubgrid> m_Grid;
};
#endif

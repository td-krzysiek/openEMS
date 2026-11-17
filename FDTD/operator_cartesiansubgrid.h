/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPERATOR_CARTESIANSUBGRID_H
#define OPERATOR_CARTESIANSUBGRID_H

#include "cartesian_subgrid.h"
#include "operator.h"
#include <memory>

class Operator_Ext_CartesianSubgrid_Conductivity;

//! Native API adapter for the experimental compact XYZ numerical core.
class Operator_CartesianSubgrid : public Operator
{
  public:
	static Operator_CartesianSubgrid *New(const std::vector<double> &box, CartesianSubgrid::Index ratio = {2, 2, 2});
	static Operator_CartesianSubgrid *NewRegions(const std::vector<std::vector<double>> &boxes,
												 const std::vector<CartesianSubgrid::Index> &ratios);
	static void ValidateRegions(const std::vector<std::vector<double>> &boxes,
								const std::vector<CartesianSubgrid::Index> &ratios);
	~Operator_CartesianSubgrid() override;
	void SetEngineType(Engine_Types::EngineType type) override;
	void SetChildLines(const CartesianSubgrid::Axes &drawingLines);
	//! Fixed execution placement, separate from physical geometry and XML.
	void SetDomainExecution(unsigned int axis, std::vector<double> drawingCuts,
		std::vector<std::unique_ptr<GPU_Backend>> backends);
	const CartesianSubgrid &ObservationGrid() const;
	void InvalidateObservations() const { m_ObservationsValid = false; }
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
	CartesianSubgrid::Point DrawingPoint(const double *drawing) const;
	Excitation *Signal() const
	{
		return m_Exc;
	}

  private:
	Operator_CartesianSubgrid(const std::vector<std::vector<double>> &boxes,
							  const std::vector<CartesianSubgrid::Index> &ratios);
	void ValidateProperties() const;
	void PrepareSources();
	void PrepareResistors();
	std::vector<std::vector<double>> m_Boxes;
	std::vector<CartesianSubgrid::Index> m_Ratios;
	CartesianSubgrid::Axes m_ChildDrawingLines;
	std::unique_ptr<CartesianSubgrid> m_Grid;
	unsigned int m_DomainAxis = 0;
	std::vector<double> m_DomainCuts;
	std::vector<std::unique_ptr<GPU_Backend>> m_DomainBackends;
	mutable bool m_ObservationsValid = true;
	//! Non-owning handle; the operator extension list retains ownership.
	Operator_Ext_CartesianSubgrid_Conductivity *m_NativeConductivity = nullptr;
};
#endif

/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#ifndef OPERATOR_CARTESIANMULTIGRID_H
#define OPERATOR_CARTESIANMULTIGRID_H

#include "multigrid_interpolation.h"
#include "operator.h"

//! Experimental planar Cartesian 2:1 multi-grid, reusing the half-step scheduler.
/*!
  X below splitX uses every second Y line. X at/above splitX uses all Y lines.
  X and Z coordinates and the timestep are shared. This first implementation
  requires an empty interface collar and PEC outer boundaries; it is not an
  embedded three-dimensional refinement box. The parent field arrays retain
  their full extent for processing, as in the cylindrical implementation.
 */
class Operator_CartesianMultiGrid : public Operator
{
  public:
	static Operator_CartesianMultiGrid *New(double splitX, unsigned int numThreads = 0);
	virtual ~Operator_CartesianMultiGrid();
	virtual void SetEngineType(Engine_Types::EngineType type);
	virtual Engine *CreateEngine();
	virtual bool SetGeometryCSX(ContinuousStructure *geo);
	virtual int CalcECOperator(DebugFlags debugFlags = None);
	virtual void SetBoundaryCondition(int *conditions);
	virtual void SetExcitationSignal(Excitation *excitation);
	virtual void SetMaterialStoreFlags(int type, bool value);
	virtual void AddExtension(Operator_Extension *extension);
	virtual bool GetYeeCoords(int component, unsigned int pos[3], double *coords, bool dualMesh) const;
	virtual void CalcStartStopLines(unsigned int &threads, std::vector<unsigned int> &start,
									std::vector<unsigned int> &stop) const;
	virtual double GetNumberCells() const;
	virtual void ShowStat() const;
	Operator *GetCoarseOperator() const { return m_Coarse; }
	unsigned int GetSplitPos() const { return m_Split; }
	const MultiGridInterpolation &GetInterpolation() const { return m_Interpolation; }
	double GetCoupledStepNormBound() const { return m_CoupledStepNormBound; }

  protected:
	explicit Operator_CartesianMultiGrid(double splitX);
	virtual void Init();
	virtual bool SetupCSXGrid(CSRectGrid *grid);
	void SetupInterpolation();
	void CheckInterface() const;
	void FillMissingMaterialStorage();
	double CoupledStepNormBound() const;
	virtual bool Calc_ECPos(int component, const unsigned int *pos, double *circuit,
							const std::vector<CSPrimitives *> &primitives) const;

	double m_SplitX;
	unsigned int m_Split;
	Operator *m_Coarse;
	MultiGridInterpolation m_Interpolation;
	double m_CoupledStepNormBound = 0;
};

#endif

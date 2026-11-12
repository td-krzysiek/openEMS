/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#ifndef CARTESIAN_MULTIGRID_COUPLING_H
#define CARTESIAN_MULTIGRID_COUPLING_H

#include "engine_subgrid.h"

//! Host field-access coupling; no assumption about basic versus SSE storage.
class CartesianMultiGridCoupling : public Engine_SubGridCoupling
{
  public:
	CartesianMultiGridCoupling(Engine *base, Engine *coarse, const MultiGridInterpolation &interpolation,
							   const unsigned int lines[3]);
	virtual void SyncVoltages();
	virtual void SyncCurrents();
	virtual void Interpolate();

  private:
	void InterpolatePlane(unsigned int x, bool currents);
	Engine *m_Base;
	Engine *m_Coarse;
	MultiGridInterpolation m_Interpolation;
	unsigned int m_Lines[3];
};

#endif

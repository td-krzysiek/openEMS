/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIANSUBGRID_ADVANCE_H
#define CARTESIANSUBGRID_ADVANCE_H

#include "operator_extension.h"
#include "engine_extension.h"
#include "FDTD/gpu_cartesiansubgrid.h"
#include <memory>

class GPU_Backend;

//! A paired extension may own the actual compact E/H recurrence.
/*!
 * This is an execution boundary, not an additional update or post-damping hook.
 * The selected pair replaces the core half sweeps; ordinary priority-ordered
 * pre/post/apply hooks still surround them. At most one provider is permitted.
 * A provider must keep complete-step observation and fail the run on an error.
 * Native material admission remains the responsibility of its operator pair.
 */
class Engine_Ext_CartesianSubgrid_Advance : public Engine_Extension
{
  public:
	explicit Engine_Ext_CartesianSubgrid_Advance(Operator_Extension *op)
		: Engine_Extension(op)
	{
	}
	virtual void UpdateElectric(double timestep, double drive) = 0;
	virtual void UpdateMagnetic(double timestep) = 0;
	//! Explicit device implementation, never a CPU fallback behind a GPU selector.
	virtual std::unique_ptr<GPU_CartesianSubgrid> CreateDevice(GPU_Backend &backend) = 0;
	//! Called only after the returned device's complete-step ReadFields.
	virtual void ReadDeviceDiagnostics(GPU_CartesianSubgrid &device) = 0;
};

//! An operator half must explicitly opt in and construct the matching provider.
class Operator_Ext_CartesianSubgrid_Advance : public Operator_Extension
{
  public:
	explicit Operator_Ext_CartesianSubgrid_Advance(Operator *op)
		: Operator_Extension(op)
	{
	}
	Engine_Ext_CartesianSubgrid_Advance *CreateEngineExtention() override = 0;
};

#endif

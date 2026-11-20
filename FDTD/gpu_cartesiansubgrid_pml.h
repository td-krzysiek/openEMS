/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef GPU_CARTESIANSUBGRID_PML_H
#define GPU_CARTESIANSUBGRID_PML_H

#include "FDTD/gpu_backend.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

//! Created by the paired engine extension from the checked operator's derivative rows.
struct CartesianSubgridPMLData
{
	struct Rows
	{
		std::vector<std::uint32_t> start{0}, source, target, group{0};
		std::vector<float> weight, low, decay, change;
	};
	Rows electric, magnetic;
	double timestep = 0;
};

//! Shared point-operation extension, with a complete-step diagnostic read only.
class GPU_CartesianSubgridPML : public GPU_Extension
{
  public:
	virtual double ReadBoundaryWorkOut() = 0;
	virtual std::size_t StorageBytes() const = 0;
	bool NeedsFlatStep(unsigned int) const override
	{
		return true;
	}
};

std::unique_ptr<GPU_CartesianSubgridPML> GPU_CreateCartesianSubgridPML(
	GPU_Backend &backend, const GPU_CartesianSubgridFields &fields, const CartesianSubgridPMLData &data);
#endif

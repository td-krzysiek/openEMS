/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef GPU_CARTESIANSUBGRID_ENERGY_H
#define GPU_CARTESIANSUBGRID_ENERGY_H
#include "gpu_cartesiansubgrid.h"

//! Optional observation-only reductions. Backend owns allocated buffers.
class GPU_CartesianSubgridEnergy
{
  public:
	virtual ~GPU_CartesianSubgridEnergy() = default;
	virtual std::array<double, 2> Read() = 0;
	virtual std::size_t StorageBytes() const = 0;
	virtual std::uint64_t DownloadBytes() const = 0;
};
std::unique_ptr<GPU_CartesianSubgridEnergy> GPU_CreateCartesianSubgridEnergy(
	GPU_Backend &backend, GPU_CartesianSubgridFields owned,
	const std::vector<CartesianSubgridEnergyCoupling> &couplings);
#endif

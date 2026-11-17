/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef GPU_CARTESIANSUBGRID_DFT_H
#define GPU_CARTESIANSUBGRID_DFT_H
#include "gpu_cartesiansubgrid.h"
#include <complex>

//! Optional owned-state accumulator. Backend owns allocations until destruction.
//! Sampling and epoch validation belong to the execution/coordinator.
class GPU_CartesianSubgridDFT
{
  public:
	virtual ~GPU_CartesianSubgridDFT() = default;
	virtual void Accumulate(const std::vector<std::complex<float>> &weights) = 0;
	virtual void Read(std::size_t frequency, std::vector<std::complex<double>> &sums) = 0;
	virtual std::size_t StorageBytes() const = 0;
	virtual std::uint64_t DownloadBytes() const = 0;
};

//! Exact device-buffer payload; overflow/invalid counts reject before allocation.
std::size_t CartesianSubgridDFTBytes(std::size_t count, std::size_t frequencies);
void CartesianSubgridDFTWeights(const std::vector<std::complex<float>> &weights, std::size_t frequencies);
std::unique_ptr<GPU_CartesianSubgridDFT> GPU_CreateCartesianSubgridDFT(
	GPU_Backend &backend, const float *field, std::size_t owned,
	const std::vector<std::uint32_t> &indices, std::size_t frequencies);
#endif

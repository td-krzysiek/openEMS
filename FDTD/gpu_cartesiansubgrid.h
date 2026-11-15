/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef GPU_CARTESIANSUBGRID_H
#define GPU_CARTESIANSUBGRID_H

#include <cstddef>
#include "subgrid_device_precision.h"
#include <cstdint>
#include <memory>
#include <vector>

class CartesianSubgrid;
class GPU_Backend;

//! Backend-neutral upload contract for the first scoped device implementation.
/*! Expanded curls are a bounded intermediate representation. They do not
 *  allocate a global fine parent lattice, but are not structured block kernels.
 *  FP32 fields/coefficients require independent numerical qualification.
 */
struct CartesianSubgridDeviceData
{
	std::vector<float> electric, magnetic, curlWeight, transposeWeight, curlLow, transposeLow, lower;
	std::vector<std::uint32_t> curlStart{0}, curlColumn, transposeStart, transposeColumn;
	std::vector<std::uint32_t> blockStart{0}, blockEdge, bandwidth, lowerStart{0};
	std::vector<std::uint32_t> sourceEdge;
	std::vector<float> sourceWeight;
	std::vector<std::uint32_t> portStart{0}, portEdge, responseStart{0}, responseEdge;
	std::vector<float> portWeight, responseWeight, rootConductance, portCurrent, portLower;
	std::vector<std::uint32_t> sourceOffStart{0}, sourceOffEdge;
	std::vector<float> sourceOffWeight, ledger{0, 0, 0, 0};
	CartesianSubgridFloat::TailInfo factorTails;
	double portTimestep = 0;
	double stableTimestep = 0;
};

//! Device execution owned by the compact engine, created by its existing backend.
/*! A backend must reject unsupported physics, never silently advance on CPU.
 *  Steps are queued in order. ReadFields waits and refreshes the CPU observation
 *  mirror once per processing batch, not between electric/magnetic sweeps.
 *  The creating GPU_Backend must outlive this object.
 */
class GPU_CartesianSubgrid
{
  public:
	virtual ~GPU_CartesianSubgrid() = default;
	virtual void Step(double timestep, double sourceAmplitude) = 0;
	virtual void ReadFields(CartesianSubgrid &grid) = 0;
	virtual std::size_t StorageBytes() const = 0;
};
#endif

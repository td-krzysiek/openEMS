/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef GPU_CARTESIANSUBGRID_H
#define GPU_CARTESIANSUBGRID_H

#include <cstddef>
#include <array>
#include "subgrid_device_precision.h"
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

class CartesianSubgrid;
class GPU_Backend;
struct Op_SubgridCurl;

//! Borrowed normalized compact buffers; not rectangular voltage/current arrays.
struct GPU_CartesianSubgridFields
{
	float *electric = nullptr, *magnetic = nullptr;
	std::size_t electricCount = 0, magneticCount = 0;
	double fixedTimestep = 0;
};

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
	CartesianSubgridFloat::TailInfo factorTails, stateTails;
	double portTimestep = 0;
	//! A source-only worker must retain source work when another worker owns ports.
	bool trackSourceWork = false;
	double stableTimestep = 0;
};

//! Paired-physics recurrence hooks around the compact algebraic solve.
//! The extension supplies the prepared factor/port data, never post-damps a step.
class GPU_CartesianSubgridRecurrence
{
  public:
	virtual ~GPU_CartesianSubgridRecurrence() = default;
	virtual void Bind(GPU_Backend &backend, const GPU_CartesianSubgridFields &fields, float *rhs) = 0;
	//! After original transpose-curl/source capture, before the prepared mass solve.
	virtual void PrepareElectric(double timestep) = 0;
	//! After the implicit/port solve, before the ordered soft-source kicks.
	virtual void FinishElectric(double timestep) = 0;
	virtual void UpdateMagnetic(const Op_SubgridCurl &curl, double timestep) = 0;
	virtual void CheckObservation() const = 0;
	virtual std::array<double, 2> ReadLedger() = 0;
	virtual std::size_t StorageBytes() const = 0;
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
	//! Paired device extensions require ordered half steps and a fixed recurrence step.
	virtual void UpdateElectric(double, double)
	{
		throw std::logic_error("CartesianSubgrid: no device electric half-step");
	}
	virtual void UpdateMagnetic(double)
	{
		throw std::logic_error("CartesianSubgrid: no device magnetic half-step");
	}
	virtual void RequireFixedTimestep(double)
	{
		throw std::logic_error("CartesianSubgrid: no fixed device timestep");
	}
	virtual GPU_CartesianSubgridFields Fields() const
	{
		return {};
	}
	//! Domain-only transport completion, after the matching ghost upload.
	virtual void AcknowledgeGhosts(bool, std::size_t)
	{
		throw std::logic_error("CartesianSubgrid: ghost acknowledgement is unavailable");
	}
	//! Complete-step loss/source-work observation; blocking like ReadFields.
	virtual std::array<double, 2> ReadLedger()
	{
		throw std::logic_error("CartesianSubgrid: device ledger observation is unavailable");
	}
	//! Extension-specific complete-step diagnostics; absent on the lossless core.
	virtual std::array<double, 2> ReadAdvanceLedger()
	{
		throw std::logic_error("CartesianSubgrid: recurrence diagnostics unavailable");
	}
};

std::unique_ptr<GPU_CartesianSubgrid> GPU_CreateCartesianSubgridAdvanced(
	GPU_Backend &backend, const CartesianSubgrid &grid, CartesianSubgridDeviceData data,
	std::unique_ptr<GPU_CartesianSubgridRecurrence> recurrence, double timestep);
#endif

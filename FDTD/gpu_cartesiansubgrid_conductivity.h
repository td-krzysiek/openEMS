/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef GPU_CARTESIANSUBGRID_CONDUCTIVITY_H
#define GPU_CARTESIANSUBGRID_CONDUCTIVITY_H
#include "gpu_cartesiansubgrid.h"

//! Setup-only SI-normalized tables. A factor preparation belongs to its operator.
struct CartesianSubgridConductivityData
{
	struct Rows
	{
		std::vector<std::uint32_t> start{0}, edge;
		std::vector<float> weight, low;
	};
	Rows conductance, trace;
	std::vector<double> magneticRate;
	double timestep = 0;
};
std::unique_ptr<GPU_CartesianSubgridRecurrence> GPU_CreateCartesianSubgridConductivity(
	CartesianSubgridConductivityData data);
#endif

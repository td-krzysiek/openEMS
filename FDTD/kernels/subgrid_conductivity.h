/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_CONDUCTIVITY_H
#define FDTD_KERNELS_SUBGRID_CONDUCTIVITY_H
#if !defined(__METAL_VERSION__)
#include "subgrid_update.h"
#endif

//! Per-trace nonnegative midpoint electric dissipation, before soft sources.
struct Op_SubgridElectricLoss
{
	CM_GLOBAL const uint32_t *start, *edge;
	CM_GLOBAL const float *weight, *low, *previous, *electric;
	CM_GLOBAL float *loss;
	float timestep;
	CM_FN void operator()(uint64_t row) const
	{
		float value = 0, correction = 0;
		for (uint32_t k = start[row]; k < start[row + 1]; ++k)
		{
			const float mid = .5f * previous[edge[k]] + .5f * electric[edge[k]];
			const float product = weight[k] * mid;
			subgrid_curl_add(value, correction, product);
			subgrid_curl_add(value, correction,
				subgrid_curl_product_error(weight[k], mid, product) + low[k] * mid);
		}
		value += correction;
		subgrid_add_work(loss, uint32_t(2 * row), timestep * value * value);
	}
};

//! Crank-Nicolson H update and its nonnegative per-face midpoint loss.
struct Op_SubgridMagneticLoss
{
	CM_GLOBAL float *magnetic, *loss;
	CM_GLOBAL const float *curl, *decay, *decayLow, *change, *changeLow, *rateStep;
	CM_FN void operator()(uint64_t face) const
	{
		const float previous = magnetic[face];
		float next = subgrid_field_increment(previous, decay[face], decayLow[face], previous);
		next = subgrid_field_increment(next, change[face], changeLow[face], curl[face]);
		magnetic[face] = next;
		const float midpoint = .5f * previous + .5f * next;
		subgrid_add_work(loss, uint32_t(2 * face), rateStep[face] * midpoint * midpoint);
	}
};

//! Hierarchical reduction only at observation boundaries, not on every step.
struct Op_SubgridLossReduce
{
	CM_GLOBAL const float *input;
	CM_GLOBAL float *output;
	uint32_t count;
	CM_FN void operator()(uint64_t group) const
	{
		const uint64_t first = 512 * group, stop = first + 512 < count ? first + 512 : count;
		float value = 0, correction = 0;
		for (uint64_t row = first; row < stop; ++row)
		{
			subgrid_curl_add(value, correction, input[2 * row]);
			subgrid_curl_add(value, correction, -input[2 * row + 1]);
		}
		output[2 * group] = value;
		output[2 * group + 1] = -correction;
	}
};
#endif

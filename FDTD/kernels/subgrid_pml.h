/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_PML_H
#define FDTD_KERNELS_SUBGRID_PML_H

#if !defined(__METAL_VERSION__)
#include "subgrid_update.h"
#endif

//! Paired compact-grid recurrence on normalized fields; one owned derivative row.
struct Op_SubgridPMLPre
{
	CM_GLOBAL const uint32_t *start;
	CM_GLOBAL const uint32_t *source;
	CM_GLOBAL const uint32_t *target;
	CM_GLOBAL const float *weight;
	CM_GLOBAL const float *low;
	CM_GLOBAL const float *decay;
	CM_GLOBAL const float *change;
	CM_GLOBAL const float *input;
	CM_GLOBAL const float *field;
	CM_GLOBAL float *history;
	CM_GLOBAL float *correction;
	CM_GLOBAL float *previous;
	float timestep, timestepLow;

	CM_FN void operator()(uint64_t row) const
	{
		float derivative = 0, error = 0;
		for (uint32_t k = start[row]; k < start[row + 1]; ++k)
		{
			const float value = input[source[k]];
			const float product = weight[k] * value;
			subgrid_curl_add(derivative, error, product);
			subgrid_curl_add(derivative, error,
				subgrid_curl_product_error(weight[k], value, product) + low[k] * value);
		}
		history[row] = decay[row] * history[row] + change[row] * (derivative + error);
		correction[row] = subgrid_field_increment(0, timestep, timestepLow, history[row]);
		previous[row] = field[target[row]];
	}
};

//! One point per distinct target, so both the natural-field read and final write are race-free.
struct Op_SubgridPMLPost
{
	CM_GLOBAL const uint32_t *group;
	CM_GLOBAL const uint32_t *target;
	CM_GLOBAL const float *correction;
	CM_GLOBAL const float *previous;
	CM_GLOBAL float *field;
	CM_GLOBAL float *work;
	CM_GLOBAL float *next;
	uint32_t apply;

	CM_FN void operator()(uint64_t point) const
	{
		const uint32_t begin = group[point], end = group[point + 1];
		if (apply)
		{
			field[target[begin]] = next[point];
			return;
		}
		const float natural = field[target[begin]];
		float finalValue = natural;
		for (uint32_t row = begin; row < end; ++row)
		{
			float before = natural;
			// Match the paired CPU Finish() reverse-prior summation order.
			for (uint32_t prior = row; prior > begin; --prior)
			{
				before += correction[prior - 1];
			}
			const float after = before + correction[row];
			const float priorChange = before - natural;
			subgrid_add_work(work, 2 * point,
				-.5f * correction[row] * (previous[row] + priorChange + after));
			finalValue = after;
		}
		next[point] = finalValue;
	}
};
#endif

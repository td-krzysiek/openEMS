/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_ENERGY_H
#define FDTD_KERNELS_SUBGRID_ENERGY_H
#if !defined(__METAL_VERSION__)
#include "subgrid_update.h"
#endif

//! Observation-only maximum, also used hierarchically on its partial outputs.
struct Op_SubgridEnergyMax
{
	CM_GLOBAL const float *input;
	CM_GLOBAL float *output;
	uint32_t count;
	CM_FN void operator()(uint64_t group) const
	{
		const uint64_t first = 512 * group, stop = first + 512 < count ? first + 512 : count;
		float maximum = 0;
		for (uint64_t row = first; row < stop; ++row)
		{
			const float x = input[row], value = x < 0 ? -x : x;
			if (value != value)
			{
				output[group] = value;
				return;
			}
			maximum = value > maximum ? value : maximum;
		}
		output[group] = maximum;
	}
};

//! Galerkin quadratic form, scaled to avoid squaring tiny/large FP32 fields.
//! Pair terms are unique off-diagonal mass entries; H has no pair entries.
struct Op_SubgridEnergyTerms
{
	CM_GLOBAL const float *input;
	CM_GLOBAL const uint32_t *first, *second;
	CM_GLOBAL const float *weight, *low;
	CM_GLOBAL float *output;
	uint32_t states, pairs;
	float scale;
	CM_FN void operator()(uint64_t group) const
	{
		const uint64_t begin = 512 * group, count = uint64_t(states) + pairs;
		const uint64_t stop = begin + 512 < count ? begin + 512 : count;
		float value = 0, correction = 0;
		for (uint64_t row = begin; row < stop; ++row)
		{
			if (row < states)
			{
				const float x = input[row] / scale, product = x * x;
				subgrid_curl_add(value, correction, .5f * product);
				subgrid_curl_add(value, correction, .5f * subgrid_curl_product_error(x, x, product));
			}
			else
			{
				const uint32_t pair = uint32_t(row - states);
				const float x = input[first[pair]] / scale, y = input[second[pair]] / scale;
				const float product = x * y, residual = subgrid_curl_product_error(x, y, product);
				const float term = weight[pair] * product;
				subgrid_curl_add(value, correction, term);
				subgrid_curl_add(value, correction, subgrid_curl_product_error(weight[pair], product, term));
				subgrid_curl_add(value, correction, low[pair] * product);
				subgrid_curl_add(value, correction, weight[pair] * residual);
				subgrid_curl_add(value, correction, low[pair] * residual);
			}
		}
		output[2 * group] = value;
		output[2 * group + 1] = -correction;
	}
};
#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_UPDATE_H
#define FDTD_KERNELS_SUBGRID_UPDATE_H

#if !defined(__METAL_VERSION__)
#include "cell_update.h"
#endif

//! Normalized curl/transpose row. The gather transpose needs no atomics.
struct Op_SubgridCurl
{
	CM_GLOBAL const uint32_t *start;
	CM_GLOBAL const uint32_t *column;
	CM_GLOBAL const float *weight;
	CM_GLOBAL const float *input;
	CM_GLOBAL float *output;
	float scale;
	uint32_t additive;

	CM_FN void operator()(uint64_t row) const
	{
		float value = 0;
		for (uint32_t k = start[row]; k < start[row + 1]; ++k)
		{
			value += weight[k] * input[column[k]];
		}
		if (additive)
		{
			output[row] += scale * value;
		}
		else
		{
			output[row] = value;
		}
	}
};

//! One independent banded interface block per point, no global field scatter races.
struct Op_SubgridMass
{
	CM_GLOBAL const uint32_t *start;
	CM_GLOBAL const uint32_t *edge;
	CM_GLOBAL const uint32_t *bandwidth;
	CM_GLOBAL const uint32_t *factorStart;
	CM_GLOBAL const float *lower;
	CM_GLOBAL float *rhs;

	CM_FN void operator()(uint64_t block) const
	{
		const uint32_t first = start[block], n = start[block + 1] - first;
		const uint32_t band = bandwidth[block], stride = band + 1, offset = factorStart[block];
		for (uint32_t i = 0; i < n; ++i)
		{
			float value = rhs[edge[first + i]];
			for (uint32_t j = i > band ? i - band : 0; j < i; ++j)
			{
				value -= lower[offset + i * stride + i - j] * rhs[edge[first + j]];
			}
			rhs[edge[first + i]] = value / lower[offset + i * stride];
		}
		for (uint32_t k = n; k > 0; --k)
		{
			const uint32_t i = k - 1;
			float value = rhs[edge[first + i]];
			const uint32_t stop = n < i + band + 1 ? n : i + band + 1;
			for (uint32_t j = i + 1; j < stop; ++j)
			{
				value -= lower[offset + j * stride + j - i] * rhs[edge[first + j]];
			}
			rhs[edge[first + i]] = value / lower[offset + i * stride];
		}
	}
};

struct Op_SubgridElectric
{
	CM_GLOBAL float *electric;
	CM_GLOBAL const float *increment;
	float timestep;

	CM_FN void operator()(uint64_t edge) const
	{
		electric[edge] += timestep * increment[edge];
	}
};

//! Dispatched serially when source indices repeat, preserving additive source order.
struct Op_SubgridSource
{
	CM_GLOBAL float *electric;
	CM_GLOBAL const uint32_t *edge;
	CM_GLOBAL const float *weight;
	float amplitude;

	CM_FN void operator()(uint64_t source) const
	{
		electric[edge[source]] += amplitude * weight[source];
	}
};
#endif

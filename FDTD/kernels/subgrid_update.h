/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_UPDATE_H
#define FDTD_KERNELS_SUBGRID_UPDATE_H

#if !defined(__METAL_VERSION__)
#include "cell_update.h"
#include <cmath>
#endif

//! Compensated sum of split-coefficient products, using only device FP32.
/*! Explicit FMA recovers product rounding; it does not enable contraction of
 *  unrelated physics. Neumaier accumulation also retains low residuals when
 *  subsequent, much larger terms cancel. Keep Metal safe-math enabled.
 */
CM_FN void subgrid_curl_add(CM_PRIVATE float &sum, CM_PRIVATE float &correction, float term)
{
	const float total = sum + term;
	const float magnitude = sum < 0 ? -sum : sum;
	const float other = term < 0 ? -term : term;
	correction += magnitude >= other ? (sum - total) + term : (term - total) + sum;
	sum = total;
}

CM_FN float subgrid_curl_product_error(float coefficient, float field, float product)
{
#if defined(__METAL_VERSION__) || defined(__HIPCC__) || defined(__CUDACC__)
	return fma(coefficient, field, -product);
#else
	return std::fma(coefficient, field, -product);
#endif
}

//! Preserve the supplied timestep and the final state addition in FP32 arithmetic.
CM_FN float subgrid_field_increment(float previous, float high, float low, float increment)
{
	const float product = high * increment;
	float value = previous, correction = 0;
	subgrid_curl_add(value, correction, product);
	subgrid_curl_add(value, correction,
		subgrid_curl_product_error(high, increment, product) + low * increment);
	return value + correction;
}

//! Normalized curl/transpose row. The gather transpose needs no atomics.
struct Op_SubgridCurl
{
	CM_GLOBAL const uint32_t *start;
	CM_GLOBAL const uint32_t *column;
	CM_GLOBAL const uint32_t *base;
	CM_GLOBAL const uint32_t *pattern;
	CM_GLOBAL const float *weight;
	CM_GLOBAL const float *low;
	CM_GLOBAL const float *input;
	CM_GLOBAL float *output;
	float scale, scaleLow;
	uint32_t additive;

	CM_FN void operator()(uint64_t row) const
	{
		const uint32_t p = pattern[row], origin = base[row];
		float value = 0, correction = 0;
		for (uint32_t k = start[p]; k < start[p + 1]; ++k)
		{
			const float field = input[origin + column[k]];
			const float product = weight[k] * field;
			const float residual = subgrid_curl_product_error(weight[k], field, product) + low[k] * field;
			subgrid_curl_add(value, correction, product);
			subgrid_curl_add(value, correction, residual);
		}
		value += correction;
		if (additive)
		{
			output[row] = subgrid_field_increment(output[row], scale, scaleLow, value);
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
	float timestep, timestepLow;

	CM_FN void operator()(uint64_t edge) const
	{
		electric[edge] = subgrid_field_increment(electric[edge], timestep, timestepLow, increment[edge]);
	}
};


CM_FN void subgrid_add_work(CM_GLOBAL float *ledger, uint32_t slot, float value)
{
	const float corrected = value - ledger[slot + 1];
	const float total = ledger[slot] + corrected;
	ledger[slot + 1] = (total - ledger[slot]) - corrected;
	ledger[slot] = total;
}

//! Preserve the unsolved transpose-curl only at source edges, before the mass solve.
struct Op_SubgridSourceCurl
{
	CM_GLOBAL const float *increment;
	CM_GLOBAL const uint32_t *edge;
	CM_GLOBAL float *curl;

	CM_FN void operator()(uint64_t source) const
	{
		curl[source] = increment[edge[source]];
	}
};

//! One small coupled Woodbury system per invocation; ports retain their CPU order.
struct Op_SubgridPortSolve
{
	CM_GLOBAL const uint32_t *start;
	CM_GLOBAL const uint32_t *edge;
	CM_GLOBAL const float *weight;
	CM_GLOBAL const float *rootConductance;
	CM_GLOBAL const float *current;
	CM_GLOBAL const float *lower;
	CM_GLOBAL const float *electric;
	CM_GLOBAL const float *increment;
	CM_GLOBAL float *rhs;
	CM_GLOBAL float *ledger;
	float timestep, amplitude;
	uint32_t count;

	CM_FN void operator()(uint64_t point) const
	{
		if (point != 0)
		{
			return;
		}
		for (uint32_t i = 0; i < count; ++i)
		{
			float voltage = 0;
			for (uint32_t k = start[i]; k < start[i + 1]; ++k)
			{
				voltage += weight[k] * (electric[edge[k]] + .5f * timestep * increment[edge[k]]);
			}
			rhs[i] = amplitude * current[i] - rootConductance[i] * voltage;
			for (uint32_t j = 0; j < i; ++j)
			{
				rhs[i] -= lower[i * count + j] * rhs[j];
			}
			rhs[i] /= lower[i * count + i];
		}
		for (uint32_t k = count; k > 0; --k)
		{
			const uint32_t i = k - 1;
			for (uint32_t j = i + 1; j < count; ++j)
			{
				rhs[i] -= lower[j * count + i] * rhs[j];
			}
			rhs[i] /= lower[i * count + i];
			const float midpoint = amplitude * current[i] - rhs[i];
			subgrid_add_work(ledger, 0, timestep * midpoint * midpoint);
			subgrid_add_work(ledger, 2, timestep * midpoint * amplitude * current[i]);
		}
	}
};

//! Dispatch serially: sparse responses of different ports may share a field edge.
struct Op_SubgridPortApply
{
	CM_GLOBAL float *electric;
	CM_GLOBAL const uint32_t *start;
	CM_GLOBAL const uint32_t *edge;
	CM_GLOBAL const float *weight;
	CM_GLOBAL const float *rhs;
	float timestep;

	CM_FN void operator()(uint64_t port) const
	{
		for (uint32_t k = start[port]; k < start[port + 1]; ++k)
		{
			electric[edge[k]] += timestep * weight[k] * rhs[port];
		}
	}
};

//! Ordered soft sources, including the exact midpoint/source-work contribution.
struct Op_SubgridSource
{
	CM_GLOBAL float *electric;
	CM_GLOBAL const uint32_t *edge;
	CM_GLOBAL const float *weight;
	CM_GLOBAL const uint32_t *offStart;
	CM_GLOBAL const uint32_t *offEdge;
	CM_GLOBAL const float *offWeight;
	CM_GLOBAL const float *sourceCurl;
	CM_GLOBAL float *ledger;
	float amplitude, timestep;
	uint32_t hasPorts;

	CM_FN void operator()(uint64_t source) const
	{
		const float change = amplitude * weight[source];
		if (hasPorts)
		{
			float derivative = electric[edge[source]];
			for (uint32_t k = offStart[source]; k < offStart[source + 1]; ++k)
			{
				derivative += offWeight[k] * electric[offEdge[k]];
			}
			subgrid_add_work(ledger, 2,
				change * (derivative + .5f * change - .5f * timestep * sourceCurl[source]));
		}
		electric[edge[source]] += change;
	}
};

#endif

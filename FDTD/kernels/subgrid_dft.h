/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_DFT_H
#define FDTD_KERNELS_SUBGRID_DFT_H
#if !defined(__METAL_VERSION__)
#include "subgrid_arithmetic.h"
#endif

//! One owned normalized state per thread; frequency-major compensated sums.
//! No ghost sampling, inter-device atomics, field downloads or device FP64.
struct Op_SubgridDFT
{
	CM_GLOBAL const float *field;
	CM_GLOBAL const uint32_t *indices;
	CM_GLOBAL float *sums;
	uint32_t count, frequency;
	float real, imag;

	CM_FN void operator()(uint64_t sample) const
	{
		const float value = field[indices[sample]];
		const uint64_t offset = 4 * (uint64_t(frequency) * count + sample);
		for (uint32_t part = 0; part < 2; ++part)
		{
			const float weight = part ? imag : real;
			const float product = weight * value;
			float high = sums[offset + 2 * part], low = sums[offset + 2 * part + 1];
			subgrid_curl_add(high, low, product);
			subgrid_curl_add(high, low, subgrid_curl_product_error(weight, value, product));
			// Keep the low part a residual after every sample. Otherwise a long
			// cancellation sequence accumulates all small terms in one FP32 low.
			const float total = high + low;
			const float recovered = total - high;
			low = (high - (total - recovered)) + (low - recovered);
			high = total;
			sums[offset + 2 * part] = high;
			sums[offset + 2 * part + 1] = low;
		}
	}
};
#endif

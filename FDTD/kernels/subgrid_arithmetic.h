/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_ARITHMETIC_H
#define FDTD_KERNELS_SUBGRID_ARITHMETIC_H
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

#endif

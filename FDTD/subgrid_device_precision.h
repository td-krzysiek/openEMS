/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DEVICE_PRECISION_H
#define CARTESIAN_SUBGRID_DEVICE_PRECISION_H

#include "subgrid_field_state.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace CartesianSubgridFloat
{
struct TailInfo
{
	std::size_t terms = 0;
	double l1 = 0, relativeBound = 0;
};

//! Scalar conversion only; the complete state still requires StateTo's bound.
inline float StateValue(double value)
{
	const float converted = static_cast<float>(value);
	if (!std::isfinite(value) || !std::isfinite(converted))
	{
		throw std::invalid_argument("CartesianSubgrid: overflowing/nonfinite initial device state");
	}
	return value != 0 && std::abs(converted) < std::numeric_limits<float>::min() ? 0 : converted;
}

//! Normal initial fields plus explicitly bounded negligible FTZ tails.
/*! Do not erase an entirely subnormal state. The summed dropped magnitude,
 *  relative to the original maximum, must be at most eps(float)^2. This makes
 *  compact Gaussian initial conditions usable without filtering resolved fields.
 */
//! Preserve global scan order and FTZ qualification without requiring storage.
template <class Sink> inline void StateTo(const std::vector<double> &values, TailInfo &info, Sink sink)
{
	double maximum = 0;
	TailInfo dropped;
	for (double value : values)
	{
		const float converted = StateValue(value);
		maximum = std::max(maximum, std::abs(value));
		if (value != 0 && std::abs(converted) < std::numeric_limits<float>::min())
		{
			++dropped.terms;
			dropped.l1 += std::abs(value);
		}
		sink(converted);
	}
	dropped.relativeBound = maximum == 0 ? 0 : dropped.l1 / maximum;
	const double epsilon = std::numeric_limits<float>::epsilon();
	if (!std::isfinite(dropped.relativeBound) || dropped.relativeBound > epsilon * epsilon)
	{
		throw std::invalid_argument("CartesianSubgrid: initial state FTZ perturbation exceeds bound");
	}
	info.terms += dropped.terms;
	info.l1 += dropped.l1;
	info.relativeBound = std::max(info.relativeBound, dropped.relativeBound);
}

inline std::vector<float> State(const std::vector<double> &values, TailInfo &info)
{
	std::vector<float> result;
	result.reserve(values.size());
	StateTo(values, info, [&](float value) { result.push_back(value); });
	return result;
}

//! Exact implicit +0 contributes no dropped FTZ tails and needs no host array.
template <class Sink> inline void StateTo(const CartesianSubgridFieldState &values, TailInfo &info, Sink sink)
{
	if (!values.ImplicitZero())
	{
		StateTo(values.Dense(), info, sink);
		return;
	}
	for (std::size_t i = 0; i < values.size(); ++i)
	{
		sink(0.f);
	}
}
inline std::vector<float> State(const CartesianSubgridFieldState &values, TailInfo &info)
{
	std::vector<float> result;
	result.reserve(values.size());
	StateTo(values, info, [&](float value) { result.push_back(value); });
	return result;
}

//! Explicit FTZ representation for negligible fill tails; never drop a diagonal.
/*! Metal/device arithmetic need not preserve FP32 subnormals. A factor may
 *  contain exponentially small band fill far below FP32 rounding. Preserve all
 *  normal values and positive diagonals; admit dropped tails only when their
 *  summed absolute perturbation is below eps(float)^2 times the smallest
 *  factor diagonal. This is a representation bound, not a physics tolerance.
 */
inline std::vector<float> Factor(const std::vector<double> &lower, std::size_t rows,
								std::size_t band, TailInfo &info)
{
	if (!rows || band >= rows || rows > std::numeric_limits<std::size_t>::max() / (band + 1) ||
		lower.size() != rows * (band + 1))
	{
		throw std::invalid_argument("CartesianSubgrid: invalid device band factor");
	}
	double minimum = std::numeric_limits<double>::infinity();
	for (std::size_t i = 0; i < rows; ++i)
	{
		const double diagonal = lower[i * (band + 1)];
		if (!(diagonal > 0) || !std::isfinite(diagonal))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid device factor diagonal");
		}
		minimum = std::min(minimum, diagonal);
	}
	std::vector<float> result;
	TailInfo dropped;
	for (std::size_t k = 0; k < lower.size(); ++k)
	{
		const double value = lower[k];
		const float converted = static_cast<float>(value);
		if (!std::isfinite(value) || !std::isfinite(converted))
		{
			throw std::invalid_argument("CartesianSubgrid: overflowing/nonfinite device factor");
		}
		if (value != 0 && std::abs(converted) < std::numeric_limits<float>::min())
		{
			if (k % (band + 1) == 0)
			{
				throw std::invalid_argument("CartesianSubgrid: device factor diagonal underflow");
			}
			++dropped.terms;
			dropped.l1 += std::abs(value);
			result.push_back(0);
		}
		else
		{
			result.push_back(converted);
		}
	}
	dropped.relativeBound = dropped.l1 / minimum;
	const double epsilon = std::numeric_limits<float>::epsilon();
	if (!std::isfinite(dropped.relativeBound) || dropped.relativeBound > epsilon * epsilon)
	{
		throw std::invalid_argument("CartesianSubgrid: device factor FTZ perturbation exceeds bound");
	}
	info.terms += dropped.terms;
	info.l1 += dropped.l1;
	info.relativeBound = std::max(info.relativeBound, dropped.relativeBound);
	return result;
}
} // namespace CartesianSubgridFloat
#endif

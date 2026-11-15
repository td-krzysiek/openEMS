/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DEVICE_PRECISION_H
#define CARTESIAN_SUBGRID_DEVICE_PRECISION_H

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

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_MATERIAL_QUADRATURE_H
#define SUBGRID_MATERIAL_QUADRATURE_H

#include <array>
#include <cmath>
#include <stdexcept>

//! Explicit scalar leaf-volume approximation. The caller resolves native CAD
//! priorities at every point. This is not a conformal dielectric-interface law.
namespace CartesianSubgridMaterialQuadrature
{
using Point = std::array<double, 3>;
//! Absolute epsilon, absolute mu, electric and magnetic SI conductivity.
using Values = std::array<double, 4>;

inline void Validate(unsigned int samples)
{
	if (samples != 0 && samples != 2 && samples != 4 && samples != 8)
	{
		throw std::invalid_argument("CartesianSubgrid: material samples must be 0, 2, 4 or 8 per axis");
	}
}

template <class Lookup> Values Average(const Point &lo, const Point &hi,
	unsigned int samples, const Lookup &lookup)
{
	Validate(samples);
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (!std::isfinite(lo[n]) || !std::isfinite(hi[n]) || !(hi[n] > lo[n]))
		{
			throw std::invalid_argument("CartesianSubgrid: finite positive material integration cell required");
		}
	}
	const unsigned int count = samples ? samples : 1;
	Values first{}, sum{}, compensation{};
	bool constant = true, initialized = false;
	for (unsigned int i = 0; i < count; ++i)
	{
		for (unsigned int j = 0; j < count; ++j)
		{
			for (unsigned int k = 0; k < count; ++k)
			{
				const unsigned int index[3] = {i, j, k};
				Point p;
				for (unsigned int n = 0; n < 3; ++n)
				{
					// Retain the previous exact centre expression for strict boxes.
					p[n] = samples ? lo[n] + (hi[n] - lo[n]) * ((index[n] + .5) / count) : (lo[n] + hi[n]) / 2;
					if (!std::isfinite(p[n]) || !(p[n] > lo[n]) || !(p[n] < hi[n]))
					{
						throw std::invalid_argument("CartesianSubgrid: material quadrature point is not representable");
					}
				}
				const Values value = lookup(p);
				for (unsigned int n = 0; n < 4; ++n)
				{
					if (!std::isfinite(value[n]) || (n < 2 ? !(value[n] > 0) : value[n] < 0))
					{
						throw std::invalid_argument("CartesianSubgrid: invalid sampled material coefficient");
					}
					// Positive dyadic weights, with compensated summation. Constant
					// leaves return the original bits, including conductivity.
					const double term = value[n] / double(count * count * count) - compensation[n];
					const double next = sum[n] + term;
					compensation[n] = (next - sum[n]) - term;
					sum[n] = next;
				}
				if (!initialized)
				{
					first = value;
					initialized = true;
				}
				constant = constant && value == first;
			}
		}
	}
	return constant ? first : sum;
}
} // namespace CartesianSubgridMaterialQuadrature
#endif

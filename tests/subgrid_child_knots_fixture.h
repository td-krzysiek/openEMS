/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_CHILD_KNOTS_FIXTURE_H
#define SUBGRID_CHILD_KNOTS_FIXTURE_H
#include "FDTD/cartesian_subgrid.h"
#include <cmath>
using Grid = CartesianSubgrid;
static constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
static Grid::Axes Base()
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int k = 0; k <= 8; ++k)
		{
			axes[n].push_back(.001 * (k + .03 * k * k + .17 * n));
		}
	}
	return axes;
}

static Grid::Axes Children(const Grid::Axes &base, const Grid::Index &ratio)
{
	Grid::Axes result;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (std::size_t k = 0; k + 1 < base[n].size(); ++k)
		{
			result[n].push_back(base[n][k]);
			for (unsigned int j = 1; j < ratio[n]; ++j)
			{
				// Deliberately unequal child lengths: equal subdivision must fail this oracle.
				const double fraction = std::pow(double(j) / ratio[n], 1.5 + .1 * n);
				result[n].push_back(base[n][k] + fraction * (base[n][k + 1] - base[n][k]));
			}
		}
		result[n].push_back(base[n].back());
	}
	return result;
}

static Grid Fixture(const Grid::Index &ratio, bool multiple, bool fine, Grid::CurlLayout layout,
	const Grid::Axes &children)
{
	const auto base = Base();
	std::vector<Grid::Region> regions{{{{2, 2, 2}, {3, 4, 4}}, ratio}};
	if (multiple)
	{
		regions.push_back({{{5, 3, 3}, {6, 6, 6}}, {1, ratio[1], ratio[2]}});
	}
	return Grid(base, epsilon, mu, regions, fine, {}, {}, {}, {}, layout, children);
}

#endif

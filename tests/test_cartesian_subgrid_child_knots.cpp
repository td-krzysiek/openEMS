/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_cartesiansubgrid.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using Grid = CartesianSubgrid;

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

#include "subgrid_child_knots_fixture.h"

static void Geometry(Grid &grid, const Grid::Axes &children)
{
	double volume = 0;
	for (std::size_t i = 0; i < grid.CellCount(); ++i)
	{
		const auto leaf = grid.LeafAddresses()[i];
		const auto lo = grid.LeafStart(i), length = grid.LeafLengths(i);
		Grid::Point middle;
		for (unsigned int n = 0; n < 3; ++n)
		{
			Require(lo[n] == children[n].at(leaf.lo[n]), "native leaf start moved");
			Require(length[n] == children[n].at(leaf.lo[n] + leaf.size[n]) - lo[n],
				"native leaf length differs from exact endpoints");
			middle[n] = lo[n] + .5 * length[n];
		}
		Require(grid.LeafAt(middle) == i, "native midpoint lookup selects another leaf");
		volume += length[0] * length[1] * length[2];
	}
	const auto base = Base();
	const double expected = (base[0].back() - base[0].front()) *
		(base[1].back() - base[1].front()) * (base[2].back() - base[2].front());
	Require(std::abs(volume / expected - 1) < 2e-13, "native children do not tile the domain volume");
	// Native child coordinates and parents are owning values, not borrowed views.
	const auto copy = grid;
	Require(copy.LeafStart(copy.CellCount() / 2) == grid.LeafStart(grid.CellCount() / 2),
		"native child coordinates were not copied");
}

static void Dynamics(const Grid::Index &ratio, bool multiple, bool fine)
{
	const auto base = Base(), children = Children(base, ratio);
	auto patterns = Fixture(ratio, multiple, fine, Grid::CurlLayout::Patterns, children);
	auto spatial = Fixture(ratio, multiple, fine, Grid::CurlLayout::SpatialBlocks, children);
	Geometry(patterns, children);
	Geometry(spatial, children);
	const auto initial = [&](unsigned int n, const Grid::Point &point)
	{
		const auto a = (n + 1) % 3, b = (n + 2) % 3;
		return std::sin(3.141592653589793 * (point[a] - base[a].front()) / (base[a].back() - base[a].front())) *
			std::sin(3.141592653589793 * (point[b] - base[b].front()) / (base[b].back() - base[b].front()));
	};
	const auto divergenceField = [](unsigned int n, const Grid::Point &p)
	{
		return std::sin(731 * p[(n + 1) % 3]) + .2 * std::cos(451 * p[(n + 2) % 3]);
	};
	Require(patterns.CurlDivergenceResidual(divergenceField) < 5e-12, "native-coordinate curl loses D C = 0");
	patterns.SetElectric(initial);
	spatial.SetElectric(initial);
	const double dt = .3 * std::min(patterns.StableTimestep(), spatial.StableTimestep());
	const double original = patterns.Energy(dt);
	double worstEnergy = 0, worstField = 0;
	for (unsigned int step = 0; step < 256; ++step)
	{
		patterns.Step(dt);
		spatial.Step(dt);
		worstEnergy = std::max(worstEnergy, std::abs(patterns.Energy(dt) / original - 1));
	}
	const auto p = patterns.DeviceStateData(), s = spatial.DeviceStateData();
	Require(p.electric.size() == s.electric.size() && p.magnetic.size() == s.magnetic.size(),
		"native structured layout changes field dimensions");
	for (std::size_t k = 0; k < p.electric.size(); ++k)
	{
		worstField = std::max(worstField, double(std::abs(p.electric[k] - s.electric[k])));
	}
	for (std::size_t k = 0; k < p.magnetic.size(); ++k)
	{
		worstField = std::max(worstField, double(std::abs(p.magnetic[k] - s.magnetic[k])));
	}
	Require(std::isfinite(worstEnergy) && worstEnergy < 2e-11, "native-coordinate leapfrog energy drift");
	Require(worstField < 1e-12, "native structured/reference fields differ");
	std::cout << "CHILD_KNOTS ratio=" << ratio[0] << '/' << ratio[1] << '/' << ratio[2]
		<< " multiple=" << multiple << " fine=" << fine << " leaves=" << patterns.CellCount()
		<< " energy_error=" << worstEnergy << " field_error=" << worstField << std::endl;
}

static void Rejections()
{
	const Grid::Index ratio{2, 2, 2};
	for (unsigned int mode = 0; mode < 5; ++mode)
	{
		auto children = Children(Base(), ratio);
		if (mode == 0)
		{
			children[0].pop_back();
		}
		if (mode == 1)
		{
			children[1][4] += 1e-9;
		}
		if (mode == 2)
		{
			children[2][3] = children[2][2];
		}
		if (mode == 3)
		{
			children[0][1] = std::numeric_limits<double>::quiet_NaN();
		}
		if (mode == 4)
		{
			children[0][1] = std::nextafter(children[0][0], children[0][2]);
		}
		bool rejected = false;
		try
		{
			Fixture(ratio, false, false, Grid::CurlLayout::Patterns, children);
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected, "invalid native child coordinates accepted");
	}
}

int main()
{
	try
	{
		for (const Grid::Index ratio : {Grid::Index{2, 1, 2}, Grid::Index{2, 4, 1}, Grid::Index{4, 2, 2}})
		{
			for (bool multiple : {false, true})
			{
				for (bool fine : {false, true})
				{
					Dynamics(ratio, multiple, fine);
				}
			}
		}
		Rejections();
		std::cout << "PASS: exact child knots, owning geometry, stable dynamics and five negative controls" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

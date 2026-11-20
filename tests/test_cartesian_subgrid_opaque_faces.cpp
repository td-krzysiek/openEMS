/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using Grid = CartesianSubgrid;
static constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
static constexpr double bottom = .0008, top = .000818;
static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}
static bool InLayer(const Grid::Point &lo, const Grid::Point &hi, unsigned int normal)
{
	const double z = .5 * (lo[normal] + hi[normal]);
	return z > bottom && z < top;
}

static void Plane(unsigned int normal, unsigned int mode, bool multiple)
{
	const unsigned int a = (normal + 1) % 3, b = (normal + 2) % 3;
	Grid::Axes lines;
	lines[normal] = {0, .0005, bottom, top, .0015, .002};
	for (unsigned int i = 0; i <= (multiple ? 8u : 4u); ++i)
	{
		lines[a].push_back(.001 * i);
	}
	lines[b] = {0, .001, .002, .003, .004};
	Grid::Region first;
	first.box = {{1, 1, 1}, {3, 3, 3}};
	first.box.stop[normal] = 4;
	first.ratio = {1, 1, 1};
	first.ratio[a] = mode == 0 ? 1 : 2;
	first.ratio[b] = mode == 0 ? 1 : 4;
	std::vector<Grid::Region> regions{first};
	if (multiple)
	{
		auto second = first;
		second.box.start[a] = 4;
		second.box.stop[a] = 6;
		std::swap(second.ratio[a], second.ratio[b]);
		regions.push_back(second);
	}
	Grid grid(lines, epsilon, mu, regions, mode == 1, {}, {}, {},
		[=](const Grid::Point &lo, const Grid::Point &hi) { return InLayer(lo, hi, normal); });
	const auto faces = grid.OpaqueBoundaryFaces();
	std::set<std::tuple<std::size_t, unsigned int, unsigned int>> unique;
	double area = 0;
	for (const auto &face : faces)
	{
		Require(face.normal == normal && grid.LeafIsOpaque(face.leaf), "internal opaque face was exposed");
		Require(unique.emplace(face.leaf, face.normal, face.side).second, "duplicate opaque boundary face");
		const auto lo = grid.LeafStart(face.leaf), length = grid.LeafLengths(face.leaf);
		const double coordinate = lo[normal] + face.side * length[normal];
		Require(std::abs(coordinate - (face.side ? top : bottom)) < 1e-17,
			"opaque surface moved away from the physical slab");
		area += length[a] * length[b];
	}
	const double expected = 2 * lines[a].back() * lines[b].back();
	Require(std::abs(area - expected) < expected * 1e-12, "opaque surface area was lost or duplicated");
	if (mode == 2)
	{
		Require(faces.size() == (multiple ? 176 : 88), "unexpected local-plane face count");
	}
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(731 * p[0] + 423 * p[1] - 511 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return std::cos(933 * p[0] - 331 * p[1] + 211 * p[2] + n) / 376.730313668;
	});
	const double dt = .5 * grid.StableTimestep(), initial = grid.Energy(dt);
	for (unsigned int step = 0; step < 256; ++step)
	{
		grid.Step(dt, 0);
	}
	const double drift = std::abs(grid.Energy(dt) - initial) / initial;
	Require(drift < 2e-9, "opaque face extraction changed the field invariant");
	std::cout << "OPAQUE_INTERNAL_FACE normal=" << normal << " mode=" << mode << " multiple=" << multiple
		<< " faces=" << faces.size() << " area=" << area << " drift=" << drift << std::endl;
}

static void PartialFace()
{
	// Deliberately make every electric trace in the thin layer PEC, so the
	// geometry reaches surface extraction instead of failing the earlier PEC
	// trace guard. A fine corner is void although the coarse-face centre is
	// opaque. Testing only that centre would silently lose the partial face.
	const Grid::Axes lines = {{{0, .001, .002, .003, .004}, {0, .001, .002, .003, .004},
		{0, .0005, bottom, top, .0015, .002}}};
	Grid::Region region;
	region.box = {{1, 1, 1}, {3, 3, 4}};
	region.ratio = {2, 2, 1};
	const auto pec = [](unsigned int, const Grid::Point &lo, const Grid::Point &hi)
	{
		return lo[2] >= bottom && hi[2] <= top;
	};
	const auto opaque = [](const Grid::Point &lo, const Grid::Point &hi)
	{
		const double x = .5 * (lo[0] + hi[0]), y = .5 * (lo[1] + hi[1]);
		const bool hole = x > .001 && x < .0015 && y > .001 && y < .0015;
		return InLayer(lo, hi, 2) && !hole;
	};
	Grid grid(lines, epsilon, mu, {region}, false, {}, pec, {}, opaque);
	try
	{
		grid.OpaqueBoundaryFaces();
	}
	catch (const std::invalid_argument &error)
	{
		Require(std::string(error.what()).find("face subdivision") != std::string::npos,
			"partial opaque face failed for an unrelated reason");
		std::cout << "PASS_REJECT_PARTIAL_OPAQUE_FACE" << std::endl;
		return;
	}
	throw std::runtime_error("partial coarse/fine exposed face was silently admitted");
}

int main()
{
	try
	{
		for (unsigned int normal = 0; normal < 3; ++normal)
		{
			for (unsigned int mode = 0; mode < 3; ++mode)
			{
				Plane(normal, mode, false);
			}
			Plane(normal, 2, true);
		}
		PartialFace();
		std::cout << "PASS: internal opaque interfaces and retained partial-face guard" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

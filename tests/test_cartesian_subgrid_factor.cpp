/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

static void Require(bool condition, const char *message)
{
	if (!condition)
		throw std::runtime_error(message);
}

static Grid Strip(unsigned int n)
{
	return Grid({n, 8, 8}, {.001, .0013, .0009}, epsilon, mu,
				{{{{2, 2, 2}, {n - 2, 6, 6}}, {2, 2, 2}}});
}

static double Seed(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

static void DenseReference()
{
	// Frozen fb98953 dense-core oracle, generated before the banded change by
	// benchmark/subgrid_core.cpp: strip 64 200, O3 C++17, double states.
	const std::array<double, 18> expected{
		-1.1167872687063207, 1.0184022343413954, .1244222426463905, .37375866028209093,
		-.45141616489152059, -.099816672104001788, .062330999771516925, .3541027536620594,
		-.82598865194840421, 2.1231991762100026, .45825823469578331, -.439297815640879,
		.18835908257520345, 1.1610901928127013, -.031095235517277588, -.4755101739450156,
		-.47578699685269787, -.84146960056945641};
	Grid grid = Strip(64);
	const auto info = grid.InterfaceStorage();
	Require(info.blocks > 1 && info.largestBlock >= 100, "long coupled interface fixture was not exercised");
	Require(info.factorEntries * 8 < info.denseEquivalentEntries, "interface factors retained dense storage");
	Require(info.workspaceEntries == info.largestBlock, "interface workspace is not block-sized");
	Require(std::abs(grid.StableTimestep() / 7.4572218243463747e-13 - 1) < 2e-13,
			"banded factor changed the conservative timestep bound");
	grid.SetElectric(Seed);
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return Seed(n, p) / 376.730313668; });
	const double dt = .8 * grid.StableTimestep(), initial = grid.Energy(dt);
	for (unsigned int step = 0; step < 200; ++step)
		grid.Step(dt);
	std::size_t index = 0;
	double error = 0;
	for (const Grid::Point &fraction : {Grid::Point{.17, .31, .43}, Grid::Point{.53, .61, .67},
									  Grid::Point{.81, .73, .29}})
	{
		const Grid::Point point{fraction[0] * 64 * .001, fraction[1] * 8 * .0013, fraction[2] * 8 * .0009};
		for (unsigned int n = 0; n < 3; ++n)
			for (double value : {grid.Electric(n, point), 376.730313668 * grid.Magnetic(n, point)})
				error = std::max(error, std::abs(value - expected[index++]));
	}
	Require(error < 5e-11, "banded solve disagrees with frozen dense E/H fields");
	for (unsigned int step = 0; step < 10000; ++step)
		grid.Step(dt);
	Require(std::abs(grid.Energy(dt) / initial - 1) < 1e-10, "banded solve broke the lossless invariant");
	std::cout << "dense-oracle maximum E/H error " << error << ", factor entries " << info.factorEntries
			  << '/' << info.denseEquivalentEntries << '\n';
}

static void AllocationGuard()
{
	// Keep the previous limit explicit: the representation improvement is not
	// qualification for arbitrarily large interfaces or a reason to raise it.
	Grid supported = Strip(512);
	Require(supported.InterfaceStorage().largestBlock == 1024, "largest supported interface was not exercised");
	Require(supported.InterfaceStorage().workspaceEntries == 1024, "workspace scales with the global field");
	bool rejected = false;
	try
	{
		Grid invalid = Strip(514);
	}
	catch (const std::invalid_argument &error)
	{
		rejected = std::string(error.what()).find("1024 edges") != std::string::npos;
	}
	Require(rejected, "oversized interface escaped the unchanged guard");
	Grid ordinary({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu, false, {1, 1, 1});
	Require(ordinary.InterfaceStorage().blocks == 0 && ordinary.InterfaceStorage().workspaceEntries == 0,
			"uncoupled grid allocated an interface solve workspace");
	ordinary.SetElectric(Seed);
	ordinary.Step(.8 * ordinary.StableTimestep());
}

int main()
{
	try
	{
		DenseReference();
		AllocationGuard();
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

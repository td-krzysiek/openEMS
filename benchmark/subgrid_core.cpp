/* SPDX-License-Identifier: GPL-3.0-or-later */
// Standalone compact-reference benchmark; no claim about ordinary Yee/GPU speed.
#include "cartesian_subgrid.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

using Grid = CartesianSubgrid;
using Clock = std::chrono::steady_clock;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

static double Seed(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

int main(int argc, char **argv)
{
	try
	{
		if (argc != 4)
			throw std::invalid_argument("usage: subgrid_core uniform|local|strip|multiple SIZE STEPS");
		const std::string kind = argv[1];
		const unsigned int n = std::stoul(argv[2]), steps = std::stoul(argv[3]);
		if (n < 8 || n > 512 || !steps)
			throw std::invalid_argument("invalid benchmark dimensions");
		Grid::Index cells{n, n, n};
		std::vector<Grid::Region> regions;
		if (kind == "strip")
		{
			cells = {n, 8, 8};
			regions = {{{{2, 2, 2}, {n - 2, 6, 6}}, {2, 2, 2}}};
		}
		else if (kind == "local" || kind == "uniform")
			regions = {{{{n / 4, n / 4, n / 4}, {3 * n / 4, 3 * n / 4, 3 * n / 4}},
						kind == "local" ? Grid::Index{2, 2, 2} : Grid::Index{1, 1, 1}}};
		else if (kind == "multiple")
			regions = {{{{1, 2, 2}, {n / 2 - 1, n - 2, n - 2}}, {2, 2, 2}},
					   {{{n / 2 + 1, 2, 2}, {n - 1, n - 2, n - 2}}, {4, 1, 2}}};
		else
			throw std::invalid_argument("unknown benchmark geometry");
		const auto start = Clock::now();
		Grid grid(cells, {.001, .0013, .0009}, epsilon, mu, regions);
		const double setup = std::chrono::duration<double>(Clock::now() - start).count();
		grid.SetElectric(Seed);
		grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return Seed(n, p) / 376.730313668; });
		const double dt = .8 * grid.StableTimestep(), initial = grid.Energy(dt);
		const auto advance = Clock::now();
		for (unsigned int i = 0; i < steps; ++i)
			grid.Step(dt);
		const double elapsed = std::chrono::duration<double>(Clock::now() - advance).count();
		double rss = 0;
#if defined(__unix__) || defined(__APPLE__)
		struct rusage usage;
		if (getrusage(RUSAGE_SELF, &usage) == 0)
#if defined(__APPLE__)
			rss = usage.ru_maxrss;
#else
			rss = usage.ru_maxrss * 1024.;
#endif
#endif
		std::cout << std::setprecision(17) << "{\"case\":\"" << kind << "\",\"size\":" << n
				  << ",\"steps\":" << steps << ",\"leaves\":" << grid.CellCount()
				  << ",\"edges\":" << grid.EdgeCount() << ",\"faces\":" << grid.FaceCount()
				  << ",\"setup_seconds\":" << setup << ",\"update_seconds\":" << elapsed
				  << ",\"accounted_bytes\":" << grid.StorageBytes() << ",\"peak_rss_bytes\":" << rss
				  << ",\"stable_timestep\":" << grid.StableTimestep()
				  << ",\"invariant_drift\":" << std::abs(grid.Energy(dt) / initial - 1);
#ifdef SUBGRID_FACTOR_STATS
		const auto info = grid.InterfaceStorage();
		std::cout << ",\"interface_blocks\":" << info.blocks << ",\"largest_block\":" << info.largestBlock
				  << ",\"largest_half_bandwidth\":" << info.largestHalfBandwidth
				  << ",\"factor_entries\":" << info.factorEntries
				  << ",\"dense_equivalent_entries\":" << info.denseEquivalentEntries
				  << ",\"interface_work_entries\":" << info.workspaceEntries;
		std::cout << ",\"inverse_bound_solves\":" << info.inverseBoundSolves
				  << ",\"bipartite_bound_blocks\":" << info.bipartiteBoundBlocks;
#endif
		std::cout << ",\"fields\":[";
		bool first = true;
		for (const Grid::Point &fraction : {Grid::Point{.17, .31, .43}, Grid::Point{.53, .61, .67},
										  Grid::Point{.81, .73, .29}})
		{
			const Grid::Point point{fraction[0] * cells[0] * .001, fraction[1] * cells[1] * .0013,
								   fraction[2] * cells[2] * .0009};
			for (unsigned int component = 0; component < 3; ++component)
				for (double value : {grid.Electric(component, point), 376.730313668 * grid.Magnetic(component, point)})
				{
					if (!first)
						std::cout << ',';
					first = false;
					std::cout << value;
				}
		}
		std::cout << "]}\n";
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

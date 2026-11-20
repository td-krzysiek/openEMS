/* SPDX-License-Identifier: GPL-3.0-or-later */
// Domain preparation only: no FDTD throughput or ordinary-Yee speedup claim.
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

int main(int argc, char **argv)
{
	try
	{
		if (argc != 4)
		{
			throw std::invalid_argument("usage: subgrid_domain_setup streamed|expanded SIZE WORKERS");
		}
		const std::string mode = argv[1];
		const unsigned int n = std::stoul(argv[2]), workers = std::stoul(argv[3]);
		if ((mode != "streamed" && mode != "expanded") || n < 8 || n > 64 || workers < 1 || workers > 8)
		{
			throw std::invalid_argument("invalid bounded setup benchmark");
		}
		using Grid = CartesianSubgrid;
		using Plan = CartesianSubgridDomainPlan;
		const auto started = std::chrono::steady_clock::now();
		Grid grid({n, n, n}, {.001, .0013, .0009}, 8.85418781762e-12, 1.2566370614359173e-6,
			std::vector<Grid::Region>{{{{n / 4, n / 4, n / 4}, {3 * n / 4, 3 * n / 4, 3 * n / 4}}, {2, 2, 2}}});
		const auto prepared = std::chrono::steady_clock::now();
		std::vector<double> cuts;
		for (unsigned int i = 1; i < workers; ++i)
		{
			cuts.push_back(grid.BaseLines()[0].back() * i / workers);
		}
		Plan plan(grid, 0, cuts, mode == "expanded" ? Plan::CurlPreparation::ExpandedOracle : Plan::CurlPreparation::Streamed);
		const auto completed = std::chrono::steady_clock::now();
		struct rusage usage;
		if (getrusage(RUSAGE_SELF, &usage) != 0)
		{
			throw std::runtime_error("peak RSS unavailable");
		}
		std::size_t rss = usage.ru_maxrss;
#ifndef __APPLE__
		rss *= 1024;
#endif
		std::size_t terms = 0, owned = 0, ghosts = 0;
		for (const auto &domain : plan.Domains())
		{
			if (domain.operators)
			{
				for (const auto *curl : {&domain.operators->curl, &domain.operators->transpose})
				{
					for (auto pattern : curl->pattern)
					{
						terms += curl->start[pattern + 1] - curl->start[pattern];
					}
				}
			}
			else
			{
				terms += domain.fields.curlColumn.size() + domain.fields.transposeColumn.size();
			}
			owned += domain.ownedElectric + domain.ownedMagnetic;
			ghosts += domain.globalElectric.size() + domain.globalMagnetic.size() - domain.ownedElectric - domain.ownedMagnetic;
		}
		std::cout << std::setprecision(17) << "{\"mode\":\"" << mode << "\",\"size\":" << n
			<< ",\"workers\":" << workers << ",\"leaves\":" << grid.CellCount()
			<< ",\"owned\":" << owned << ",\"ghosts\":" << ghosts << ",\"local_terms\":" << terms
			<< ",\"grid_bytes\":" << grid.StorageBytes() << ",\"peak_rss_bytes\":" << rss
			<< ",\"grid_seconds\":" << std::chrono::duration<double>(prepared - started).count()
			<< ",\"plan_seconds\":" << std::chrono::duration<double>(completed - prepared).count()
			<< ",\"timestep\":" << grid.StableTimestep() << "}" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

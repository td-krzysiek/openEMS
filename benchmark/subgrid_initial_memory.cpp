/* SPDX-License-Identifier: GPL-3.0-or-later */
// Count/host-memory characterization, not an end-to-end GPU benchmark.
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

using Clock = std::chrono::steady_clock;
using Grid = CartesianSubgrid;

static std::size_t PeakRSS()
{
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage) != 0)
	{
		throw std::runtime_error("process peak RSS unavailable");
	}
#ifdef __APPLE__
	return usage.ru_maxrss;
#else
	return std::size_t(usage.ru_maxrss) * 1024;
#endif
}

int main(int argc, char **argv)
{
	try
	{
		if (argc != 3)
		{
			throw std::invalid_argument("usage: characterize patterns|blocks 128|160|200");
		}
		const std::string mode(argv[1]);
		const unsigned int n = std::stoul(argv[2]);
		if ((mode != "patterns" && mode != "blocks") || (n != 128 && n != 160 && n != 200))
		{
			throw std::invalid_argument("case outside bounded characterization plan");
		}
		const std::size_t expected = std::size_t(n) * n * n + 7 * std::size_t(n / 2) * (n / 2) * (n / 2);
		Grid::Axes axes;
		const double spacing[3] = {.001, .0013, .0009};
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (unsigned int k = 0; k <= n; ++k)
			{
				axes[axis].push_back(spacing[axis] * k);
			}
		}
		const std::vector<Grid::Region> regions{{{{n / 4, n / 4, n / 4}, {3 * n / 4, 3 * n / 4, 3 * n / 4}}, {2, 2, 2}}};
		bool rejectedByDefault = false;
		if (expected > 1000000)
		{
			try
			{
				Grid rejected(axes, 8.85418781762e-12, 1.2566370614359173e-6, regions);
			}
			catch (const std::invalid_argument &error)
			{
				rejectedByDefault = std::string(error.what()).find("allocation limit exceeded") != std::string::npos;
			}
			if (!rejectedByDefault)
			{
				throw std::runtime_error("default one-million guard did not reject the larger model");
			}
		}
		const auto begin = Clock::now();
		Grid grid(axes, 8.85418781762e-12, 1.2566370614359173e-6, regions, false, {}, {}, {}, {},
			mode == "blocks" ? Grid::CurlLayout::SpatialBlocks : Grid::CurlLayout::Patterns, {}, 16000000);
		if (grid.CellCount() != expected)
		{
			throw std::runtime_error("larger topology count differs from independent region arithmetic");
		}

		const auto built = Clock::now();
		const auto gridRSS = PeakRSS();
		const auto initialBytes = grid.StorageBytes();
		std::size_t owned = 0, ghosts = 0, packages = 0;
		const std::vector<double> cuts{axes[0][n / 4], axes[0][n / 2], axes[0][3 * n / 4]};
		CartesianSubgridDomainPlan plan(grid, 0, cuts, CartesianSubgridDomainPlan::CurlPreparation::Streamed,
			[&](std::size_t worker, const CartesianSubgridDomainData &domain)
			{
				if (worker != packages || !domain.operators)
				{
					throw std::runtime_error("consumed worker package order/operator mismatch");
				}
				++packages;
				owned += domain.ownedElectric + domain.ownedMagnetic;
				ghosts += domain.globalElectric.size() + domain.globalMagnetic.size() - domain.ownedElectric - domain.ownedMagnetic;
			});
		const auto prepared = Clock::now();
		const auto planRSS = PeakRSS();
		if (packages != 4 || plan.UploadsRetained() || owned != grid.EdgeCount() + grid.FaceCount())
		{
			throw std::runtime_error("larger ownership or consumed-payload identity failed");
		}
		const double dt = .7 * grid.StableTimestep(), initial = grid.Energy(dt);
		grid.Step(dt);
		const double drift = std::abs(grid.Energy(dt));
		if (initial != 0)
		{
			throw std::runtime_error("zero-state memory control was preseeded");
		}
		if (!std::isfinite(drift) || drift > 1e-10)
		{
			throw std::runtime_error("larger first-step energy invariant failed");
		}
		const auto completed = Clock::now();
		const auto storage = grid.StructuredCurlStorage();
		std::cout << std::setprecision(17) << "{\"mode\":\"" << mode << "\",\"size\":" << n
			<< ",\"leaves\":" << grid.CellCount() << ",\"expected_leaves\":" << expected
			<< ",\"default_rejected\":" << (rejectedByDefault ? "true" : "false")
			<< ",\"leaf_budget\":16000000,\"edges\":" << grid.EdgeCount() << ",\"faces\":" << grid.FaceCount()
			<< ",\"owned\":" << owned << ",\"ghosts\":" << ghosts << ",\"packages\":" << packages
			<< ",\"grid_bytes_initial\":" << initialBytes
			<< ",\"grid_bytes_after_step\":" << grid.StorageBytes() << ",\"grid_peak_rss\":" << gridRSS
			<< ",\"plan_peak_rss\":" << planRSS << ",\"after_step_peak_rss\":" << PeakRSS()
			<< ",\"grid_seconds\":" << std::chrono::duration<double>(built - begin).count()
			<< ",\"plan_seconds\":" << std::chrono::duration<double>(prepared - built).count()
			<< ",\"warm_step_and_energy_seconds\":" << std::chrono::duration<double>(completed - prepared).count()
			<< ",\"timestep\":" << dt << ",\"first_step_energy_drift\":" << drift
			<< ",\"regular_rows\":" << storage.regularRows << ",\"exception_rows\":" << storage.exceptionRows
			<< ",\"field_samples\":[";
		for (unsigned int k = 0; k < 12; ++k)
		{
			const Grid::Point p{axes[0].back() * (.13 + .051 * k), axes[1].back() * (.21 + .037 * k),
				axes[2].back() * (.17 + .043 * k)};
			if (k)
			{
				std::cout << ',';
			}
			std::cout << '[' << grid.Electric(k % 3, p) << ',' << grid.Magnetic(k % 3, p) << ']';
		}
		std::cout << "],\"curl_samples\":[";
		std::vector<std::pair<std::size_t, double>> row;
		for (std::size_t k = 0; k < 64; ++k)
		{
			const auto face = (grid.FaceCount() - 1) * k / 63;
			grid.CurlRow(face, row);
			if (k)
			{
				std::cout << ',';
			}
			std::cout << '[' << face;
			for (const auto &term : row)
			{
				std::cout << ",[" << term.first << ',' << term.second << ']';
			}
			std::cout << ']';
		}
		std::cout << "]}" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

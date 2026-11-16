/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using Grid = CartesianSubgrid;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

template <class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const std::invalid_argument &)
	{
		return;
	}
	throw std::runtime_error("invalid reference sweep accepted");
}

static void Control(bool ports)
{
	Grid grid({12, 10, 8}, {.001, .0013, .0009}, 8.85418781762e-12, 1.2566370614359173e-6,
		std::vector<Grid::Region>{{{{3, 3, 2}, {9, 7, 6}}, {2, 2, 2}}});
	const double dt = .8 * grid.StableTimestep();
	if (ports)
	{
		grid.SetResistors({{{.001, .0013, .0009}, {.002, .0013, .0009}, 50, .02}}, dt);
	}
	grid.SetSources({{0, .02}});
	grid.SetElectric([](unsigned int n, const Grid::Point &point)
	{
		return std::sin(1700 * point[0] + 900 * point[1] - 2100 * point[2] + n);
	});
	const auto cold = grid.StorageBytes();
	{
		// Neither state export nor cold whole-domain preparation uses CPU work.
		const auto state = grid.DeviceStateData();
		const auto curls = grid.DeviceCurls(false);
		CartesianSubgridDomainPlan plan(grid, 0, {.006});
		const auto again = grid.DeviceStateData();
		Require(state.electric == again.electric && state.magnetic == again.magnetic,
			"device preparation changed host field state");
		Require(!curls.curl.pattern.empty() && plan.Domains().size() == 2,
			"device preparation control did not prepare operators/domains");
		Require(grid.StorageBytes() == cold, "device preparation allocated reference sweep work");
	}
	Reject([&] { grid.UpdateElectric(2 * grid.StableTimestep()); });
	Reject([&] { grid.UpdateMagnetic(dt); });
	Require(grid.StorageBytes() == cold, "rejected sweeps allocated reference work");
	grid.UpdateElectric(dt, .01);
	const auto electric = grid.StorageBytes();
	const std::size_t electricWork = grid.EdgeCount() * sizeof(double) * (ports ? 2 : 1);
	Require(electric >= cold + electricWork, "first CPU E sweep did not allocate its workspace");
	Reject([&] { grid.UpdateElectric(dt); });
	Require(grid.StorageBytes() == electric, "rejected second E sweep changed workspace");
	grid.UpdateMagnetic(dt);
	const auto warm = grid.StorageBytes();
	Require(warm >= electric + grid.FaceCount() * sizeof(double),
		"first CPU H sweep did not allocate its workspace");
	for (unsigned int step = 0; step < 100; ++step)
	{
		grid.Step(dt, .001 * std::sin(.1 * step));
		Require(grid.StorageBytes() == warm, "CPU sweep workspace grew after first use");
	}
	Require(std::isfinite(grid.Energy(dt)) && grid.Energy(dt) > 0,
		"reference workspace control lost its field state");
	std::cout << "REFERENCE_WORK ports=" << ports << " cold_bytes=" << cold
		<< " first_cpu_work_bytes=" << warm - cold << " repeat_capacity_stable=1" << std::endl;
}

int main()
{
	try
	{
		Control(false);
		Control(true);
		std::cout << "PASS: device preparation needs no CPU sweep scratch; reference sweeps allocate once" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

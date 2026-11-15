/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_cartesiansubgrid.h"
#include "subgrid_device_target.h"
#include "subgrid_fp32_oracle.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

static double Seed(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

static Grid Fixture(unsigned int axis)
{
	Grid::Index cells{8, 8, 8}, stop{6, 6, 6};
	cells[axis] = 768;
	stop[axis] = 766;
	Grid grid(cells, {.001, .0013, .0009}, epsilon, mu,
			  {{{{2, 2, 2}, stop}, {2, 2, 2}}});
	Require(grid.InterfaceStorage().largestBlock == 1536, "large interface fixture did not cross old limit");
	grid.SetElectric(Seed);
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return Seed(n, p) / 376.730313668; });
	return grid;
}

static double Error(const std::vector<float> &a, const std::vector<float> &b)
{
	Require(a.size() == b.size(), "large device field size differs");
	double difference = 0, scale = 0;
	for (std::size_t i = 0; i < a.size(); ++i)
	{
		Require(std::isfinite(b[i]), "nonfinite large device field");
		difference = std::max(difference, std::abs(double(a[i]) - b[i]));
		scale = std::max(scale, std::abs(double(a[i])));
	}
	return difference / std::max(scale, 1e-300);
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool driven)
{
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "selected backend unavailable");
	target.Verify(*backend);
	auto cpu = Fixture(axis), physical = Fixture(axis), oracle = Fixture(axis);
	if (driven)
	{
		const std::size_t edge = cpu.EdgeCount() / 2;
		for (auto *grid : {&cpu, &physical, &oracle})
		{
			grid->SetSources({{edge, .1}, {edge + 1, -.03}});
		}
	}
	auto execution = backend->CreateCartesianSubgrid(physical);
	auto reference = CreateCartesianSubgridReference(oracle);
	Require(bool(execution) && bool(reference), "large device execution unavailable");
	const double dt = .8 * cpu.StableTimestep(), initial = cpu.Energy(dt);
	for (unsigned int step = 0; step < 512; ++step)
	{
		const double drive = driven ? std::exp(-std::pow((double(step) - 20) / 8, 2)) : 0;
		cpu.Step(dt, drive);
		execution->Step(dt, drive);
		reference->Step(dt, drive);
		if ((step + 1) % 128 == 0)
		{
			execution->ReadFields(physical);
			reference->ReadFields(oracle);
		}
	}
	const auto a = cpu.DeviceData(), b = physical.DeviceData(), c = oracle.DeviceData();
	const double error = std::max(Error(a.electric, b.electric), Error(a.magnetic, b.magnetic));
	const double oracleError = std::max(Error(c.electric, b.electric), Error(c.magnetic, b.magnetic));
	const double drift = std::abs(physical.Energy(dt) / initial - 1);
	Require(error < 2e-4, "large device/double field error exceeded 2e-4");
	Require(oracleError < 2e-6, "large device/FP32 oracle error exceeded 2e-6");
	if (!driven)
	{
		Require(drift < 5e-4, "large device invariant drift exceeded 5e-4");
	}
	std::cout << "axis=" << axis << " driven=" << driven << " CPU_error=" << error
			  << " FP32_oracle_error=" << oracleError << " invariant_drift=" << drift
			  << " device_bytes=" << execution->StorageBytes() << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "usage: test_cartesian_subgrid_device_large [--metal|--cuda]");
		const SubgridDeviceTarget target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, false);
		}
		Control(target, 0, true);
		std::cout << "PASS: bounded large interfaces on selected backend" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

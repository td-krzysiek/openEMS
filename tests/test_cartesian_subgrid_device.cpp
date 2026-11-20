/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_cartesiansubgrid.h"
#include "subgrid_fp32_oracle.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

static double Seed(unsigned int n, const Grid::Point &point)
{
	return std::sin(1731 * point[0] + 933 * point[1] - 2671 * point[2] + n) +
		   .3 * std::cos(4381 * point[0] - 2271 * point[1] + 411 * point[2] + 2 * n);
}

static Grid Fixture(const std::string &kind)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		axes[n].push_back(0);
		for (unsigned int i = 0; i < 12; ++i)
		{
			axes[n].push_back(axes[n].back() + .001 * (1 + .1 * n + (kind == "graded" ? .02 * i : 0)));
		}
	}
	Grid::Index ratio = kind == "uniform" ? Grid::Index{1, 1, 1} :
		kind == "graded" ? Grid::Index{2, 4, 2} : Grid::Index{2, 2, 2};
	Grid grid(axes, epsilon, mu, {{{{3, 3, 3}, {9, 9, 9}}, ratio}});
	grid.SetElectric(Seed);
	grid.SetMagnetic([](unsigned int n, const Grid::Point &point) { return Seed(n, point) / 376.730313668; });
	return grid;
}

static double Error(const std::vector<float> &reference, const std::vector<float> &candidate)
{
	Require(reference.size() == candidate.size(), "device field size differs");
	double delta = 0, scale = 0;
	for (std::size_t i = 0; i < reference.size(); ++i)
	{
		Require(std::isfinite(candidate[i]), "nonfinite device field");
		delta = std::max(delta, std::abs(double(candidate[i]) - reference[i]));
		scale = std::max(scale, std::abs(double(reference[i])));
	}
	return delta / std::max(scale, 1e-300);
}

static void Control(const SubgridDeviceTarget &target, const std::string &kind, unsigned int steps, bool source)
{
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "selected backend unavailable during repeated allocation");
	target.Verify(*backend);
	auto cpu = Fixture(kind), device = Fixture(kind), oracle = Fixture(kind);
	if (source)
	{
		const std::size_t edge = cpu.EdgeCount() / 2;
		for (auto *grid : {&cpu, &device, &oracle})
		{
			grid->SetSources({{edge, .1}, {edge, -.03}, {edge + 1, .08}});
		}
	}
	auto execution = backend->CreateCartesianSubgrid(device);
	auto reference = CreateCartesianSubgridReference(oracle);
	Require(bool(execution) && bool(reference), "selected backend does not implement subgridding");
	Require(execution->StorageBytes() > 0, "no device allocation reported");
	const double dt = .8 * cpu.StableTimestep(), initial = cpu.Energy(dt);
	for (unsigned int step = 0; step < steps; ++step)
	{
		const double drive = source ? std::exp(-std::pow((double(step) - 20) / 8, 2)) : 0;
		cpu.Step(dt, drive);
		execution->Step(dt, drive);
		reference->Step(dt, drive);
		if ((step + 1) % 128 == 0)
		{
			// Same processing boundary as a native batch: no field upload/restart.
			execution->ReadFields(device);
			reference->ReadFields(oracle);
		}
	}
	execution->ReadFields(device);
	reference->ReadFields(oracle);
	const auto a = cpu.DeviceData(), b = device.DeviceData(), c = oracle.DeviceData();
	const double error = std::max(Error(a.electric, b.electric), Error(a.magnetic, b.magnetic));
	const double oracleError = std::max(Error(c.electric, b.electric), Error(c.magnetic, b.magnetic));
	const double drift = std::abs(device.Energy(dt) - initial) / std::abs(initial);
	std::cout << kind << " steps=" << steps << " source=" << source << " CPU_error=" << error
			  << " FP32_oracle_error=" << oracleError << " invariant_drift=" << drift
			  << " device_bytes=" << execution->StorageBytes() << std::endl;
	Require(error < 2e-4, "device/double reference field error exceeded 2e-4");
	Require(oracleError < 2e-6, "device/FP32 reference error exceeded 2e-6");
	if (!source)
	{
		Require(drift < 5e-4, "device invariant drift exceeded 5e-4");
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
	throw std::runtime_error("invalid device state was accepted");
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "usage: test_cartesian_subgrid_device [--metal|--cuda]");
		const SubgridDeviceTarget target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
		Require(bool(backend), "no selected physical backend; no fallback allowed");
		target.Verify(*backend);
		std::cout << "Selected backend: " << backend->GetName() << std::endl;
		backend.reset();
		for (const std::string kind : {"uniform", "local", "graded"})
		{
			Control(target, kind, 2000, true);
		}
		Control(target, "local", 20000, false);
		auto grid = Fixture("local");
		auto data = grid.DeviceData();
		Reject([&] { grid.ImportDeviceFields(data.electric.data(), data.electric.size() - 1,
											 data.magnetic.data(), data.magnetic.size()); });
		data.electric[0] = std::numeric_limits<float>::quiet_NaN();
		Reject([&] { grid.ImportDeviceFields(data.electric.data(), data.electric.size(),
											 data.magnetic.data(), data.magnetic.size()); });
		grid.UpdateElectric(.8 * grid.StableTimestep());
		Reject([&] { grid.DeviceData(); });
		std::cout << "PASS: scoped device fields, source ordering, FP32 oracle and complete-step guards" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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

static Grid Fixture(bool drive, bool soft)
{
	Grid grid({12, 12, 12}, {.001, .001, .001}, {{3, 3, 3}, {9, 9, 9}}, epsilon, mu);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return .1 * std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .1 * std::cos(1381 * p[0] - 2171 * p[1] + 311 * p[2] + n) / 376.730313668;
	});
	const auto initial = grid.DeviceData();
	const double dt = .8 * grid.StableTimestep();
	grid.SetResistors({{{.001, .004, .004}, {.001, .004, .005}, 50, drive ? .001 : 0},
					   {{.006, .006, .006}, {.006, .006, .0065}, 75, 0}}, dt);
	if (soft)
	{
		Require(!initial.blockEdge.empty(), "no coupled interface source edge");
		const auto edge = initial.blockEdge.front();
		// Repeated kicks on a mass-coupled edge exercise ordered work/scatter.
		grid.SetSources({{edge, 10}, {edge, -3}});
		const auto data = grid.DeviceData();
		Require(!data.sourceOffEdge.empty(), "soft source lacks mass off-diagonal neighbours");
	}
	return grid;
}

static double FieldError(const Grid &reference, const Grid &candidate)
{
	const auto a = reference.DeviceData(), b = candidate.DeviceData();
	double error = 0;
	for (const auto &fields : {std::make_pair(&a.electric, &b.electric), std::make_pair(&a.magnetic, &b.magnetic)})
	{
		double scale = 0, delta = 0;
		Require(fields.first->size() == fields.second->size(), "field dimensions differ");
		for (std::size_t i = 0; i < fields.first->size(); ++i)
		{
			scale = std::max(scale, std::abs(double(fields.first->at(i))));
			delta = std::max(delta, std::abs(double(fields.first->at(i)) - fields.second->at(i)));
		}
		error = std::max(error, delta / std::max(scale, 1e-300));
	}
	return error;
}

static void Control(const std::string &backendName, bool drive, bool soft, unsigned int steps, bool dcStress = false)
{
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(backendName));
	Require(bool(backend), "requested backend unavailable");
	if (backendName == "auto")
	{
		Require(backend->GetName().find("Metal (") == 0, "physical test requires Metal; no fallback");
	}
	auto cpu = Fixture(drive, soft), device = Fixture(drive, soft);
	auto execution = backend->CreateCartesianSubgrid(device);
	const double dt = .8 * cpu.StableTimestep(), initial = cpu.Energy(dt);
	for (unsigned int i = 0; i < steps; ++i)
	{
		const double position = (double(i) - 40) / 12;
		// Preserve the original strong-DC stress as a separate unchanged control.
		const double signal = std::exp(-position * position) * (soft && !dcStress ? position : 1);
		cpu.Step(dt, signal);
		execution->Step(dt, signal);
		if ((i + 1) % 128 == 0)
		{
			execution->ReadFields(device);
		}
	}
	execution->ReadFields(device);
	const double energyScale = std::max({std::abs(initial), std::abs(cpu.SourceWork()), cpu.ResistorLoss()});
	const double fieldError = FieldError(cpu, device);
	const double lossError = std::abs(cpu.ResistorLoss() - device.ResistorLoss()) / energyScale;
	const double workError = std::abs(cpu.SourceWork() - device.SourceWork()) / energyScale;
	const double balance = std::abs(device.Energy(dt) + device.ResistorLoss() - device.SourceWork() - initial) / energyScale;
	std::cout << backend->GetName() << " drive=" << drive << " soft=" << soft << " steps=" << steps
			  << " field_error=" << fieldError << " loss_error=" << lossError << " work_error=" << workError
			  << " balance=" << balance << std::endl;
	Require(fieldError < 2e-4, "device port fields exceed 2e-4");
	Require(lossError < 2e-4 && workError < 2e-4 && balance < 2e-4, "device port work/loss/balance exceeds 2e-4");
	Require(device.ResistorLoss() > 0, "passive resistance did not dissipate energy");
	if (!drive && !soft)
	{
		Require(device.Energy(dt) <= initial * (1 + 2e-4), "unforced passive energy grew");
	}
	try
	{
		execution->Step(dt * .5, 0);
	}
	catch (const std::invalid_argument &)
	{
		return;
	}
	throw std::runtime_error("changed fixed port timestep was accepted");
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		bool metal = false, dcStress = false;
		for (int i = 1; i < argc; ++i)
		{
			const std::string argument(argv[i]);
			if (argument == "--metal")
			{
				metal = true;
			}
			else if (argument == "--dc-stress")
			{
				dcStress = true;
			}
			else
			{
				throw std::invalid_argument("usage: test_cartesian_subgrid_device_ports [--metal] [--dc-stress]");
			}
		}
		const std::string backend = metal ? "auto" : "reference";
		if (dcStress)
		{
			Control(backend, false, true, 10000, true);
			return 0;
		}
		Control(backend, false, false, 10000);
		Control(backend, true, false, 10000);
		Control(backend, false, true, 10000);
		Control(backend, true, true, 20000);
		std::cout << "PASS: device midpoint ports, coupled soft work, passive decay and fixed-step guard" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "subgrid_fp32_oracle.h"
#include "subgrid_child_knots_fixture.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

static double Error(const std::vector<float> &expected, const std::vector<float> &actual)
{
	Require(expected.size() == actual.size(), "child-coordinate device field dimensions differ");
	double delta = 0, scale = 0;
	for (std::size_t k = 0; k < expected.size(); ++k)
	{
		Require(std::isfinite(actual[k]) && std::isfinite(expected[k]), "nonfinite child-coordinate field");
		delta = std::max(delta, std::abs(double(expected[k]) - actual[k]));
		scale = std::max(scale, std::abs(double(expected[k])));
	}
	return delta / std::max(scale, 1e-300);
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool multiple, bool ports)
{
	const auto base = Base();
	const Grid::Index ratio{2, 4, 2};
	auto grid = Fixture(ratio, multiple, false, Grid::CurlLayout::SpatialBlocks, Children(base, ratio));
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(731 * p[(n + 1) % 3]) + .2 * std::cos(451 * p[(n + 2) % 3]);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .001 * std::cos(317 * p[n]);
	});
	const double dt = .5 * grid.StableTimestep();
	if (ports)
	{
		grid.SetResistors({{{base[0][1], base[1][5], base[2][3]},
			{base[0][1], base[1][5], base[2][4]}, 50, .001}}, dt);
	}
	grid.SetSources({{grid.EdgeCount() / 3, .03}});
	auto cpu = grid, single = grid, distributed = grid, independent = grid;
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "requested child-coordinate backend unavailable");
	target.Verify(*backend);
	std::cout << "CHILD_KNOTS_SINGLE " << backend->GetName() << std::endl;
	auto one = backend->CreateCartesianSubgrid(single);
	// This independent FP32 oracle has no resistor implementation. Port
	// cases instead require CPU accuracy and exact single/domain agreement.
	auto oracle = ports ? std::unique_ptr<GPU_CartesianSubgrid>{} : CreateCartesianSubgridReference(independent);
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(grid, axis,
		std::vector<double>{base[axis][4]});
	std::vector<std::unique_ptr<GPU_Backend>> workers;
	std::vector<GPU_Backend *> pointers;
	for (std::size_t k = 0; k < plan->Domains().size(); ++k)
	{
		target.SelectDomain(k);
		workers.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(workers.back()), "requested child-coordinate domain backend unavailable");
		target.Verify(*workers.back());
		std::cout << "CHILD_KNOTS_WORKER " << k << ' ' << workers.back()->GetName() << std::endl;
		pointers.push_back(workers.back().get());
	}
	auto domains = CreateCartesianSubgridDomains(plan, pointers);
	for (unsigned int step = 0; step < 256; ++step)
	{
		const double drive = std::sin(.17 * step);
		cpu.Step(dt, drive);
		one->Step(dt, drive);
		if (oracle)
		{
			oracle->Step(dt, drive);
		}
		domains->Step(dt, drive);
	}
	one->ReadFields(single);
	if (oracle)
	{
		oracle->ReadFields(independent);
	}
	domains->ReadFields(distributed);
	const auto c = cpu.DeviceStateData(), s = single.DeviceStateData();
	const auto d = distributed.DeviceStateData(), o = independent.DeviceStateData();
	const double cpuError = std::max(Error(c.electric, s.electric), Error(c.magnetic, s.magnetic));
	const double oracleError = oracle ? std::max(Error(o.electric, s.electric), Error(o.magnetic, s.magnetic)) : -1;
	Require(cpuError < 2e-4, "child-coordinate device/CPU field error exceeds 2e-4");
	if (oracle)
	{
		Require(oracleError < 2e-6, "child-coordinate device/independent FP32 field error exceeds 2e-6");
	}
	Require(s.electric == d.electric && s.magnetic == d.magnetic,
		"child-coordinate domains change one-device FP32 fields");
	std::cout << "CHILD_KNOTS_DEVICE axis=" << axis << " multiple=" << multiple
		<< " ports=" << ports
		<< " cpu_error=" << cpuError << " oracle_error=" << oracleError
		<< " exact_domains=1" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected one optional physical backend flag");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool multiple : {false, true})
			{
				for (bool ports : {false, true})
				{
					Control(target, axis, multiple, ports);
				}
			}
		}
		std::cout << "PASS: native child coordinates on single/domain backends and independent FP32 oracle" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

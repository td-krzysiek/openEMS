/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

using Grid = CartesianSubgrid;
using Plan = CartesianSubgridDomainPlan;

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

static bool SameCurl(const CartesianSubgridDeviceCurl &a, const CartesianSubgridDeviceCurl &b)
{
	return std::tie(a.start, a.offset, a.base, a.pattern, a.weight, a.residual) ==
		std::tie(b.start, b.offset, b.base, b.pattern, b.weight, b.residual);
}

static void SamePackage(const CartesianSubgridDomainData &a, const CartesianSubgridDomainData &b)
{
	Require(a.globalElectric == b.globalElectric && a.globalMagnetic == b.globalMagnetic &&
		a.sendElectric == b.sendElectric && a.sendMagnetic == b.sendMagnetic &&
		a.ownedElectric == b.ownedElectric && a.ownedMagnetic == b.ownedMagnetic,
		"streaming changed final ownership/exports");
	const auto &x = a.fields;
	const auto &y = b.fields;
	Require(std::tie(x.electric, x.magnetic, x.lower, x.blockStart, x.blockEdge, x.bandwidth,
		x.lowerStart, x.sourceEdge, x.sourceWeight, x.portStart, x.portEdge, x.responseStart,
		x.responseEdge, x.portWeight, x.responseWeight, x.rootConductance, x.portCurrent,
		x.portLower, x.sourceOffStart, x.sourceOffEdge, x.sourceOffWeight, x.ledger,
		x.portTimestep, x.trackSourceWork, x.stableTimestep) ==
		std::tie(y.electric, y.magnetic, y.lower, y.blockStart, y.blockEdge, y.bandwidth,
		y.lowerStart, y.sourceEdge, y.sourceWeight, y.portStart, y.portEdge, y.responseStart,
		y.responseEdge, y.portWeight, y.responseWeight, y.rootConductance, y.portCurrent,
		y.portLower, y.sourceOffStart, y.sourceOffEdge, y.sourceOffWeight, y.ledger,
		y.portTimestep, y.trackSourceWork, y.stableTimestep),
		"streaming changed numeric state, factors, sources or ports");
	Require(a.operators && b.operators && SameCurl(a.operators->curl, b.operators->curl) &&
		SameCurl(a.operators->transpose, b.operators->transpose), "streaming changed exact local curls");
}

static Grid Fixture(bool ports)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int k = 0; k <= 8; ++k)
		{
			axes[n].push_back(.001 * k * (1 + .03 * n + .01 * k));
		}
	}
	Grid grid(axes, 8.85418781762e-12, 1.2566370614359173e-6,
		std::vector<Grid::Region>{{{{2, 2, 2}, {6, 6, 6}}, {2, 4, 2}}}, false);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(791 * p[0] - 1183 * p[1] + 913 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .001 * std::cos(517 * p[0] + 1233 * p[1] - 723 * p[2] + n);
	});
	if (ports)
	{
		grid.SetResistors({{{axes[0][1], axes[1][1], axes[2][1]},
			{axes[0][1], axes[1][1], axes[2][2]}, 50, 1e-5}}, .7 * grid.StableTimestep());
	}
	else
	{
		grid.SetSources({{grid.EdgeCount() / 2, .1}});
	}
	return grid;
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool ports)
{
	auto original = Fixture(ports), streamed = original;
	const double dt = .7 * original.StableTimestep();
	const auto &lines = original.BaseLines()[axis];
	const std::vector<double> cuts{lines[3], lines[5]};
	auto retained = std::make_shared<Plan>(original, axis, cuts);
	Require(retained->UploadsRetained(), "ordinary plan unexpectedly consumed");
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend *> devices;
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers;
	const CartesianSubgridDomainData *previous = nullptr;
	const auto consume = [&](std::size_t k, const CartesianSubgridDomainData &domain)
	{
		Require(k == workers.size(), "upload callback order changed");
		if (previous)
		{
			Require(previous->fields.electric.empty() && previous->fields.magnetic.empty() &&
				!previous->operators, "previous numeric package survived until the next callback");
		}
		SamePackage(domain, retained->Domains().at(k));
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "requested upload backend unavailable");
		target.Verify(*backends.back());
		devices.push_back(backends.back().get());
		workers.push_back(CreateCartesianSubgridWorker(*backends.back(), domain));
		backends.back()->Activate();
		backends.back()->Synchronize();
		previous = &domain;
	};
	auto consumed = std::make_shared<Plan>(streamed, axis, cuts, Plan::CurlPreparation::Streamed, consume);
	Require(!consumed->UploadsRetained() && workers.size() == retained->Domains().size(),
		"streamed upload completion was not recorded");
	for (const auto &domain : consumed->Domains())
	{
		Require(domain.fields.electric.empty() && domain.fields.magnetic.empty() &&
			domain.fields.lower.empty() && !domain.operators, "numeric payload retained after upload");
	}
	bool rejected = false;
	try
	{
		CreateCartesianSubgridDomains(consumed, devices);
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "a consumed plan was silently reuploaded");
	auto transport = CreateCartesianSubgridHostTransport(consumed, workers);
	auto execution = CreateCartesianSubgridCoordinator(consumed, std::move(workers), std::move(transport));
	std::vector<std::unique_ptr<GPU_Backend>> controls;
	std::vector<GPU_Backend *> controlDevices;
	for (std::size_t k = 0; k < retained->Domains().size(); ++k)
	{
		target.SelectDomain(k);
		controls.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(controls.back()), "retained control backend unavailable");
		target.Verify(*controls.back());
		controlDevices.push_back(controls.back().get());
	}
	auto baseline = CreateCartesianSubgridDomains(retained, controlDevices);
	baseline->RequireFixedTimestep(dt);
	execution->RequireFixedTimestep(dt);
	for (unsigned int step = 0; step < 128; ++step)
	{
		const double drive = std::exp(-std::pow((double(step) - 20) / 8, 2));
		baseline->Step(dt, drive);
		execution->Step(dt, drive);
	}
	baseline->ReadFields(original);
	execution->ReadFields(streamed);
	const auto a = original.DeviceStateData(), b = streamed.DeviceStateData();
	Require(a.electric == b.electric && a.magnetic == b.magnetic &&
		baseline->ReadLedger() == execution->ReadLedger(), "streamed uploads changed deterministic domain fields or work");
	Require(retained->UploadsRetained() && !retained->Domains().front().fields.electric.empty(),
		"streaming invalidated a separately retained reusable plan");
	std::cout << "UPLOAD_CONSUMER axis=" << axis << " ports=" << ports
		<< " domains=" << consumed->Domains().size() << " exact_fields=1 released_payloads=1" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "unexpected upload consumer arguments");
		const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, false);
			Control(target, axis, true);
		}
		auto grid = Fixture(false);
		std::size_t calls = 0;
		bool rejected = false;
		try
		{
			Plan aborted(grid, 0, {.003, .006}, Plan::CurlPreparation::Streamed,
				[&](std::size_t, const CartesianSubgridDomainData &)
				{
					++calls;
					throw std::invalid_argument("intentional consumer failure");
				});
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected && calls == 1, "failed upload consumer continued preparation");
		std::cout << "PASS: streamed upload lifetimes, exact fields, final exports and failure guards" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "subgrid_fp32_oracle.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
using Plan = CartesianSubgridDomainPlan;

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
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid selective observation operation accepted");
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

static std::vector<std::uint32_t> Selection(std::size_t count, unsigned int samples)
{
	std::vector<std::uint32_t> result;
	for (unsigned int k = 0; k < samples; ++k)
	{
		result.push_back(static_cast<std::uint32_t>((count - 1) * (samples - 1 - k) / (samples - 1)));
	}
	result.push_back(result.front());
	result.push_back(result[3]);
	return result;
}

static double Error(const std::vector<float> &a, const std::vector<float> &b)
{
	Require(a.size() == b.size(), "selective observation count differs");
	double delta = 0, scale = 0;
	for (std::size_t i = 0; i < a.size(); ++i)
	{
		Require(std::isfinite(a[i]) && std::isfinite(b[i]), "nonfinite selective observation");
		delta = std::max(delta, std::abs(double(a[i]) - b[i]));
		scale = std::max(scale, std::abs(double(a[i])));
	}
	return delta / std::max(scale, 1e-300);
}

static std::vector<float> Selected(const std::vector<float> &all, const std::vector<std::uint32_t> &indices)
{
	std::vector<float> result;
	for (auto index : indices)
	{
		result.push_back(all.at(index));
	}
	return result;
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool ports)
{
	auto original = Fixture(ports), one = original, split = original, oracleGrid = original;
	const double dt = .7 * original.StableTimestep();
	const auto e = Selection(original.EdgeCount(), 17), h = Selection(original.FaceCount(), 23);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> singleBackend(GPU_Backend::New(target.backend));
	Require(bool(singleBackend), "single observation backend missing");
	target.Verify(*singleBackend);
	auto single = singleBackend->CreateCartesianSubgrid(one);
	// The independent FP32 helper has no port solver. Port controls use
	// the independent double core plus the retained single/domain FP32 gate.
	auto oracle = ports ? std::unique_ptr<GPU_CartesianSubgrid>{} : CreateCartesianSubgridReference(oracleGrid);
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers;
	const auto &lines = original.BaseLines()[axis];
	const auto consume = [&](std::size_t k, const CartesianSubgridDomainData &domain)
	{
		Require(k == workers.size(), "selective worker upload order changed");
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "selective domain backend missing");
		target.Verify(*backends.back());
		std::cout << "OBSERVATION_BACKEND owner=" << k << " backend=" << backends.back()->GetName() << '\n';
		workers.push_back(CreateCartesianSubgridWorker(*backends.back(), domain));
		backends.back()->Activate();
		backends.back()->Synchronize();
	};
	auto plan = std::make_shared<Plan>(split, axis, std::vector<double>{lines[3], lines[5]},
		Plan::CurlPreparation::Streamed, consume);
	Require(!plan->UploadsRetained(), "selective test retained numeric upload packages");
	auto transport = CreateCartesianSubgridHostTransport(plan, workers);
	auto distributed = CreateCartesianSubgridCoordinator(plan, std::move(workers), std::move(transport));
	Require(single->HostObservationStorageBytes() == 0 && distributed->HostObservationStorageBytes() == 0,
		"unrequested complete host mirrors were allocated at device/coordinator startup");
	std::vector<float> ev, hv;
	for (auto *execution : {single.get(), distributed.get()})
	{
		Reject([&] { execution->ReadObservations(ev, hv); });
		Reject([&] { execution->ConfigureObservations({static_cast<std::uint32_t>(original.EdgeCount())}, h); });
		Reject([&] { execution->ConfigureObservations(e, {static_cast<std::uint32_t>(original.FaceCount())}); });
		execution->ConfigureObservations(e, h);
		Reject([&] { execution->ConfigureObservations(e, h); });
		Reject([&] { execution->ReadObservations(ev, ev); });
		execution->ReadObservations(ev, hv);
		const auto initial = original.DeviceStateData();
		Require(ev == Selected(initial.electric, e) && hv == Selected(initial.magnetic, h),
			"startup selection order, duplicates or initial values changed");
	}
	const auto singleBytes = single->HostObservationStorageBytes();
	const auto splitBytes = distributed->HostObservationStorageBytes();
	Require(singleBytes == 0 && splitBytes < 4096,
		"sparse selection materialized complete host observation storage");
	double maximum = 0, singleDomain = 0;
	for (unsigned int step = 0; step < 64; ++step)
	{
		const double x = (double(step) - 20) / 8, drive = std::exp(-x * x);
		for (auto *execution : {single.get(), distributed.get()})
		{
			execution->UpdateElectric(dt, drive);
			Reject([&] { execution->ReadObservations(ev, hv); });
			execution->UpdateMagnetic(dt);
		}
		if (ports)
		{
			oracleGrid.Step(dt, drive);
		}
		else
		{
			oracle->Step(dt, drive);
		}
		if (step % 8 == 7)
		{
			if (oracle)
			{
				oracle->ReadFields(oracleGrid);
			}
			const auto expected = oracleGrid.DeviceStateData();
			single->ReadObservations(ev, hv);
			const auto singleE = ev, singleH = hv;
			maximum = std::max({maximum, Error(Selected(expected.electric, e), ev),
				Error(Selected(expected.magnetic, h), hv)});
			distributed->ReadObservations(ev, hv);
			maximum = std::max({maximum, Error(Selected(expected.electric, e), ev),
				Error(Selected(expected.magnetic, h), hv)});
			singleDomain = std::max({singleDomain, Error(singleE, ev), Error(singleH, hv)});
			Require(single->HostObservationStorageBytes() == singleBytes &&
				distributed->HostObservationStorageBytes() == splitBytes,
				"repeated sparse reads grew hidden complete mirrors");
		}
	}
	Require(maximum < (ports ? 2e-4 : 2e-6) && singleDomain < 2e-6,
		"selective observations exceeded retained double/FP32 reference or single/domain bounds");
	for (auto pair : {std::make_pair(single.get(), &one), std::make_pair(distributed.get(), &split)})
	{
		pair.first->ReadObservations(ev, hv);
		pair.first->ReadFields(*pair.second);
		const auto complete = pair.second->DeviceStateData();
		Require(ev == Selected(complete.electric, e) && hv == Selected(complete.magnetic, h),
			"selective gather differs from exact same-epoch complete download");
		Require(pair.first->HostObservationStorageBytes() >=
			(original.EdgeCount() + original.FaceCount()) * sizeof(float),
			"explicit complete download did not account for its legacy mirror");
	}
	std::cout << "SELECTIVE_OBSERVATIONS axis=" << axis << " ports=" << ports
		<< " steps=64 samples_e=" << e.size() << " samples_h=" << h.size()
		<< " reference_error=" << maximum << " reference_bound=" << (ports ? 2e-4 : 2e-6)
		<< " single_domain_error=" << singleDomain << " single_host_bytes=" << singleBytes
		<< " domains_host_bytes=" << splitBytes << " exact_full_read=1 consumed_uploads=1\n";
}

static void Empty(const SubgridDeviceTarget &target)
{
	auto grid = Fixture(false);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "empty-selection backend missing");
	target.Verify(*backend);
	auto execution = backend->CreateCartesianSubgrid(grid);
	execution->ConfigureObservations({}, {});
	std::vector<float> e{1}, h{2};
	execution->ReadObservations(e, h);
	Require(e.empty() && h.empty() && execution->HostObservationStorageBytes() == 0,
		"empty selection read a full field or retained old outputs");
	execution->Step(.7 * grid.StableTimestep(), 0);
	execution->ReadObservations(e, h);
	Reject([&] { execution->ConfigureObservations({}, {}); });
	Require(e.empty() && h.empty() && execution->HostObservationStorageBytes() == 0,
		"empty selection changed after stepping");
	std::cout << "SELECTIVE_EMPTY steps=1 host_bytes=0\n";
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "unexpected selective-observation arguments");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		Empty(target);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, false);
			Control(target, axis, true);
		}
		std::cout << "PASS: ordered sparse normalized states, complete epochs, consumed domains and lazy full mirrors\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

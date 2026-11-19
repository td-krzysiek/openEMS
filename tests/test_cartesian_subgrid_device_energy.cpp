/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_energy_test.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
static constexpr double Epsilon = 8.85418781762e-12, Mu = 1.2566370614359173e-6;
static void Require(bool ok, const char *message)
{
	if (!ok)
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
	throw std::runtime_error("invalid device energy operation accepted");
}

static Grid Fixture(bool opaque)
{
	Grid::Axes axes, child;
	const unsigned int ratio[3] = {2, 4, 2};
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int k = 0; k <= 8; ++k)
		{
			axes[n].push_back(.001 * k * (1 + .03 * n + .01 * k));
		}
		for (unsigned int k = 0; k < 8; ++k)
		{
			for (unsigned int j = 0; j < ratio[n]; ++j)
			{
				const double fraction = j == 0 ? 0 : j == 1 ? .3 : j == 2 ? .55 : .8;
				child[n].push_back(axes[n][k] + (axes[n][k + 1] - axes[n][k]) * fraction);
			}
		}
		child[n].push_back(axes[n].back());
	}
	const auto medium = [&](const Grid::Point &lo, const Grid::Point &)
	{
		return Grid::Medium{Epsilon * (lo[0] >= axes[0][4] ? 3.25 : 1), 1.3 * Mu};
	};
	const auto metal = [&](const Grid::Point &lo, const Grid::Point &)
	{
		return opaque && lo[0] >= axes[0][7] && lo[1] >= axes[1][7] && lo[2] >= axes[2][7];
	};
	Grid grid(axes, Epsilon, Mu, std::vector<Grid::Region>{{{{2, 2, 2}, {6, 6, 6}}, {2, 4, 2}}},
		false, medium, {}, {}, metal, Grid::CurlLayout::SpatialBlocks, child);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(791 * p[0] - 1183 * p[1] + 913 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .001 * std::cos(517 * p[0] + 1233 * p[1] - 723 * p[2] + n);
	});
	return grid;
}


static void Control(const SubgridDeviceTarget &target, unsigned int axis, double scale, bool ports)
{
	auto original = Fixture(false);
	original.SetElectric([&](unsigned int n, const Grid::Point &p)
	{
		return scale * std::sin(791 * p[0] - 1183 * p[1] + 913 * p[2] + n);
	});
	original.SetMagnetic([&](unsigned int n, const Grid::Point &p)
	{
		return scale * .001 * std::cos(517 * p[0] + 1233 * p[1] - 723 * p[2] + n);
	});
	const auto &lines = original.BaseLines();
	const double dt = .7 * original.StableTimestep();
	if (ports)
	{
		original.SetResistors({{{lines[0][1], lines[1][1], lines[2][1]},
			{lines[0][1], lines[1][1], lines[2][2]}, 50, 0}}, dt);
	}
	auto one = original, split = original, unchanged = original, cpu = original;
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "energy backend missing");
	target.Verify(*backend);
	auto single = backend->CreateCartesianSubgrid(one);
	auto control = backend->CreateCartesianSubgrid(unchanged);
	Reject([&] { single->ReadEnergy(); });
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend *> pointers;
	for (unsigned int k = 0; k < 3; ++k)
	{
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "energy domain backend missing");
		target.Verify(*backends.back());
		std::cout << "ENERGY_BACKEND worker=" << k << " backend=" << backends.back()->GetName() << '\n';
		pointers.push_back(backends.back().get());
	}
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(split, axis,
		std::vector<double>{lines[axis][3], lines[axis][5]});
	auto distributed = CreateCartesianSubgridDomains(plan, pointers);
	Reject([&] { distributed->ReadEnergy(); });
	single->ConfigureEnergy(original);
	distributed->ConfigureEnergy(original);
	Reject([&] { single->ConfigureEnergy(original); });
	Reject([&] { distributed->ConfigureEnergy(original); });
	Require(single->HostObservationStorageBytes() == 0 && distributed->HostObservationStorageBytes() == 0,
		"energy configuration allocated complete host mirrors");
	const double initial = original.ElectricEnergy() + original.MagneticEnergy();
	double maximum = 0, domainError = 0;
	for (unsigned int step = 0; step < 128; ++step)
	{
		for (auto *execution : {single.get(), distributed.get()})
		{
			execution->UpdateElectric(dt, 0);
			const auto bytes = execution->EnergyDownloadBytes();
			Reject([&] { execution->ReadEnergy(); });
			Require(execution->EnergyDownloadBytes() == bytes, "half-step energy downloaded data");
			execution->UpdateMagnetic(dt);
		}
		control->Step(dt, 0);
		cpu.Step(dt, 0);
		if (step % 16 == 15)
		{
			const auto a = single->ReadEnergy(), b = distributed->ReadEnergy();
			if (step == 15)
			{
				Require(single->HostObservationStorageBytes() == 0 && distributed->HostObservationStorageBytes() == 0,
					"energy read materialized complete host fields");
			}
			const double expected[2] = {cpu.ElectricEnergy(), cpu.MagneticEnergy()};
			for (unsigned int field = 0; field < 2; ++field)
			{
				maximum = std::max(maximum, std::abs(a[field] - expected[field]) / initial);
				domainError = std::max(domainError, std::abs(a[field] - b[field]) / initial);
			}
			single->ReadFields(one);
			distributed->ReadFields(split);
			control->ReadFields(unchanged);
			CheckSubgridEnergy(*single, one, 1, "single");
			CheckSubgridEnergy(*distributed, split, 3, "domains");
			const auto x = one.DeviceStateData(), y = unchanged.DeviceStateData();
			Require(x.electric == y.electric && x.magnetic == y.magnetic,
				"energy observation changed time-stepped fields");
		}
	}
	Require(maximum < 2e-4 && domainError < 2e-6,
		"energy independent-double or single/domain gate failed");
	std::cout << "DEVICE_ENERGY axis=" << axis << " scale=" << scale << " ports=" << ports
		<< " steps=128 double_error=" << maximum << " single_domain_error=" << domainError
		<< " exact_unchanged_fields=1 pairs=" << original.ElectricEnergyCouplings().size() << '\n';
}

static void Guards(const SubgridDeviceTarget &target)
{
	auto grid = Fixture(false);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "energy guard backend missing");
	target.Verify(*backend);
	const auto count = static_cast<std::uint32_t>(grid.EdgeCount());
	const std::vector<std::vector<CartesianSubgridEnergyCoupling>> invalid{
		{{0, 0, .1}}, {{0, count, .1}}, {{0, 1, std::numeric_limits<double>::quiet_NaN()}},
		{{1, 2, .1}, {0, 1, .1}}, {{0, 1, .1}, {0, 1, .1}}, {{0, 1, 1}}};
	for (const auto &terms : invalid)
	{
		auto device = backend->CreateCartesianSubgrid(grid);
		Reject([&] { device->ConfigureEnergyCouplings(terms); });
		Reject([&] { device->ReadEnergy(); });
	}
	grid.SetElectric([](unsigned int, const Grid::Point &) { return 1.; });
	const double full = grid.ElectricEnergy();
	auto diagonal = backend->CreateCartesianSubgrid(grid);
	diagonal->ConfigureEnergyCouplings({}); // Deliberately wrong for a coupled grid.
	const double wrong = diagonal->ReadEnergy()[0];
	Require(std::abs(wrong - full) / full > 1e-4,
		"energy control cannot distinguish a diagonal-only shortcut");
	std::cout << "DEVICE_ENERGY_DIAGONAL_NEGATIVE relative_error=" << std::abs(wrong - full) / full << '\n';
	grid.SetElectric([](unsigned int, const Grid::Point &) { return 0.; });
	grid.SetMagnetic([](unsigned int, const Grid::Point &) { return 0.; });
	auto zero = backend->CreateCartesianSubgrid(grid);
	zero->ConfigureEnergy(grid);
	Require(zero->ReadEnergy() == std::array<double, 2>{0, 0}, "zero field energy differs");
	CheckSubgridEnergy(*zero, grid, 1, "zero");
	for (float value : {std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()})
	{
		auto broken = backend->CreateCartesianSubgrid(grid);
		broken->ConfigureEnergy(grid);
		const auto fields = broken->Fields();
		backend->Activate();
		backend->DeviceUpload(fields.electric, &value, sizeof(value));
		Reject([&] { broken->ReadEnergy(); });
	}
	std::cout << "DEVICE_ENERGY_GUARDS invalid_couplings=6 nonfinite_fields=2 zero=1\n";
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "unexpected device energy arguments");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		Guards(target);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, 1e-18, false);
			Control(target, axis, 1e30, false);
			Control(target, axis, 1, true);
		}
		std::cout << "PASS: scaled Galerkin field-energy reductions, domains and read-only observations\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

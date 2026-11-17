/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_observations.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
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
	throw std::runtime_error("invalid physical observation accepted");
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
		return Grid::Medium{Epsilon * (lo[0] >= axes[0][4] ? 3.556 : 1), 1.3 * Mu};
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
	grid.SetSources({{grid.EdgeCount() / 2, .1}});
	return grid;
}

static void Analytic()
{
	Grid::Axes axes;
	for (auto &a : axes)
	{
		a = {0, .001, .002, .003, .004};
	}
	Grid grid(axes, Epsilon, Mu, std::vector<Grid::Region>{{{{1, 1, 1}, {3, 3, 3}}, {2, 2, 2}}});
	grid.SetElectric([](unsigned int n, const Grid::Point &) { return 1. + n; });
	grid.SetMagnetic([](unsigned int n, const Grid::Point &) { return .001 * (1 + n); });
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (const auto &p : std::vector<Grid::Point>{{.001, .001, .001},
			{.0013, .0022, .0027}, {.003, .002, .001}})
		{
			Require(std::abs(grid.Evaluate(grid.ElectricSample(n, p)) - (1 + n)) < 3e-13,
				"analytic SI electric sample has wrong normalization");
			Require(std::abs(grid.Evaluate(grid.MagneticSample(n, p)) - .001 * (1 + n)) < 3e-16,
				"analytic SI magnetic sample has wrong normalization");
		}
	}
	Reject([&] { grid.ElectricSample(3, {.001, .001, .001}); });
	Reject([&] { grid.MagneticSample(3, {.001, .001, .001}); });
	Reject([&] { grid.ElectricSample(0, {-.001, 0, 0}); });
	Reject([&] { grid.MagneticSample(0, {0, 0, std::numeric_limits<double>::quiet_NaN()}); });
	Grid::Integral bad;
	bad.weights.push_back({grid.EdgeCount(), 1});
	Reject([&] { CartesianSubgridObservations reader(grid, {bad}); });
	bad.weights = {{0, std::numeric_limits<double>::infinity()}};
	Reject([&] { CartesianSubgridObservations reader(grid, {bad}); });
	std::cout << "PHYSICAL_ANALYTIC samples=18 invalid=6\n";
}

static double Error(const std::vector<double> &a, const std::vector<double> &b,
	const std::vector<bool> &magnetic)
{
	Require(a.size() == b.size() && a.size() == magnetic.size(), "physical observation counts differ");
	double delta[2]{}, scale[2]{};
	for (std::size_t i = 0; i < a.size(); ++i)
	{
		Require(std::isfinite(a[i]) && std::isfinite(b[i]), "nonfinite physical observation");
		const unsigned int n = magnetic[i];
		delta[n] = std::max(delta[n], std::abs(a[i] - b[i]));
		scale[n] = std::max(scale[n], std::abs(a[i]));
	}
	return std::max(delta[0] / std::max(scale[0], 1e-300), delta[1] / std::max(scale[1], 1e-300));
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool opaque)
{
	auto original = Fixture(opaque), one = original, split = original, cpu = original;
	std::vector<Grid::Point> points;
	const auto &lines = original.BaseLines();
	for (unsigned int k = 0; k <= 8; ++k)
	{
		points.push_back({lines[0][k], lines[1][k], lines[2][k]});
		if (k < 8)
		{
			Grid::Point p;
			for (unsigned int n = 0; n < 3; ++n)
			{
				p[n] = lines[n][k] + .37 * (lines[n][k + 1] - lines[n][k]);
			}
			points.push_back(p);
		}
	}
	std::vector<Grid::Integral> functionals;
	std::vector<bool> magnetic;
	for (const auto &p : points)
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			const auto e = original.ElectricSample(n, p), h = original.MagneticSample(n, p);
			Require(std::abs(original.Evaluate(e) - original.Electric(n, p)) < 3e-13 &&
				std::abs(original.Evaluate(h) - original.Magnetic(n, p)) < 3e-16,
				"physical sample differs from existing direct reconstruction");
			functionals.push_back(e);
			functionals.push_back(h);
		}
	}
	functionals.push_back(original.VoltagePath({lines[0][1], lines[1][1], lines[2][1]},
		{lines[0][1], lines[1][1], lines[2][2]}));
	functionals.push_back(original.CurrentContour(2, {lines[0][2], lines[1][2], lines[2][2]},
		{lines[0][5], lines[1][5], lines[2][2]}));
	functionals.push_back(functionals.front()); // Caller order and duplicate outputs remain.
	for (const auto &f : functionals)
	{
		magnetic.push_back(f.magnetic);
	}
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "physical observation backend missing");
	target.Verify(*backend);
	auto single = backend->CreateCartesianSubgrid(one);
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend *> borrowed;
	for (unsigned int k = 0; k < 3; ++k)
	{
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "physical observation domain backend missing");
		target.Verify(*backends.back());
		std::cout << "PHYSICAL_OBSERVER_BACKEND worker=" << k << " backend=" << backends.back()->GetName() << '\n';
		borrowed.push_back(backends.back().get());
	}
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(split, axis,
		std::vector<double>{lines[axis][3], lines[axis][5]});
	auto distributed = CreateCartesianSubgridDomains(plan, borrowed);
	CartesianSubgridObservations a(original, functionals), b(original, functionals);
	std::vector<double> av, bv, expected;
	Reject([&] { a.Read(av); });
	a.Bind(*single);
	b.Bind(*distributed);
	Reject([&] { a.Bind(*single); });
	const auto singleBytes = single->HostObservationStorageBytes();
	Require(singleBytes == 0 && distributed->HostObservationStorageBytes() < 16384,
		"physical reader requested complete mirrors");
	const double dt = .7 * original.StableTimestep();
	double reconstruction = 0, reference = 0, domains = 0;
	for (unsigned int step = 0; step < 64; ++step)
	{
		const double x = (double(step) - 20) / 8, drive = std::exp(-x * x);
		for (auto *execution : {single.get(), distributed.get()})
		{
			execution->UpdateElectric(dt, drive);
			Reject([&] { (execution == single.get() ? a : b).Read(av); });
			execution->UpdateMagnetic(dt);
		}
		cpu.Step(dt, drive);
		if (step % 8 == 7)
		{
			a.Read(av);
			b.Read(bv);
			domains = std::max(domains, Error(av, bv, magnetic));
			expected.clear();
			for (const auto &f : functionals)
			{
				expected.push_back(cpu.Evaluate(f));
			}
			reference = std::max(reference, Error(expected, av, magnetic));
			for (const auto &pair : {std::make_pair(single.get(), &one), std::make_pair(distributed.get(), &split)})
			{
				pair.first->ReadFields(*pair.second); // Qualification oracle only.
				expected.clear();
				for (const auto &f : functionals)
				{
					expected.push_back(pair.second->Evaluate(f));
				}
				reconstruction = std::max(reconstruction, Error(expected,
					pair.first == single.get() ? av : bv, magnetic));
			}
		}
	}
	Require(reconstruction < 3e-13 && reference < 2e-4 && domains < 2e-6,
		"physical E/H/U/I reconstruction or independent reference bound exceeded");
	std::cout << "PHYSICAL_OBSERVATIONS axis=" << axis << " opaque=" << opaque
		<< " steps=64 samples=" << functionals.size() << " e_states=" << a.SelectedElectricCount()
		<< " h_states=" << a.SelectedMagneticCount() << " reconstruction_error=" << reconstruction
		<< " double_error=" << reference << " single_domain_error=" << domains << '\n';
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "unexpected physical observation arguments");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		Analytic();
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, false);
			Control(target, axis, true);
		}
		std::cout << "PASS: sparse physical E/H/U/I samples and unchanged direct reconstruction\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

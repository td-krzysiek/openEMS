/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_child_knots_fixture.h"
#include "subgrid_device_target.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/subgrid_observations.h"
#include "FDTD/subgrid_device_precision.h"
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

template <class Error, class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const Error &)
	{
		return;
	}
	throw std::runtime_error("invalid implicit-zero state operation accepted");
}

static void EmptyFields(const Grid &grid)
{
	const auto info = grid.HostFieldStorage();
	Require(info.electricValues == grid.EdgeCount() && info.magneticValues == grid.FaceCount(),
		"implicit fields lost their logical sizes");
	Require(info.electricImplicitZero && info.magneticImplicitZero &&
		info.electricBytes == 0 && info.magneticBytes == 0,
		"read-only setup or sparse observation materialized full double fields");
}

static void StorageControls()
{
	CartesianSubgridFieldState huge;
	huge.AssignZero(std::size_t(1) << 30);
	const auto &read = huge;
	Require(read.size() == (std::size_t(1) << 30) && read.at(read.size() - 1) == 0 &&
		!std::signbit(read.at(0)) && read.capacity() == 0, "logical zero state allocated or changed values");
	Reject<std::logic_error>([&] { read.Dense(); });
	Reject<std::out_of_range>([&] { read.at(read.size()); });
	Reject<std::out_of_range>([&] { huge.at(huge.size()); });
	auto copied = huge;
	auto moved = std::move(copied);
	Require(moved.size() == read.size() && moved.capacity() == 0 && copied.size() == 0,
		"implicit copy/move allocated or changed size");
	CartesianSubgridFieldState small;
	small.AssignZero(7);
	Require(small.MutableDense() == std::vector<double>(7, 0.), "CPU materialization changed initial zeros");
	small.at(1) = -0.;
	small[3] = 1.25;
	small[4] = 1e-100;
	auto dense = small;
	Require(std::memcmp(dense.Dense().data(), small.Dense().data(), 7 * sizeof(double)) == 0,
		"dense copy changed exact values");
	dense[3] = 2.;
	Require(small.Dense()[3] == 1.25 && std::signbit(small.Dense()[1]), "dense state copy aliases its parent");
	CartesianSubgridFloat::TailInfo a, b;
	const auto expected = CartesianSubgridFloat::State(small.Dense(), a);
	const auto actual = CartesianSubgridFloat::State(small, b);
	Require(expected == actual && a.terms == b.terms && a.l1 == b.l1 && a.relativeBound == b.relativeBound,
		"implicit-state wrapper changed FTZ conversion");
	small.AssignZero(7);
	Require(small.capacity() == 0 && small.ImplicitZero(), "zero reset retained a dense allocation");
	CartesianSubgridFloat::TailInfo tails;
	Require(CartesianSubgridFloat::State(small, tails) == std::vector<float>(7, 0.f) && tails.terms == 0,
		"implicit zero float conversion failed");
	small[0] = 1e-50;
	Reject<std::invalid_argument>([&] { CartesianSubgridFloat::State(small, tails); });
	std::cout << "IMPLICIT_ZERO_STORAGE logical_values=1073741824 allocated_bytes=0 dense_copy=1 ftz_guards=1\n";
}

static Grid ZeroFixture(Grid::CurlLayout layout)
{
	const auto base = Base();
	const Grid::Index ratio{2, 4, 2};
	return Fixture(ratio, true, false, layout, Children(base, ratio));
}

static void DenseZeros(Grid &grid)
{
	grid.SetElectric([](unsigned int, const Grid::Point &) { return 0.; });
	grid.SetMagnetic([](unsigned int, const Grid::Point &) { return 0.; });
}

static void CPUControl(Grid::CurlLayout layout)
{
	auto lazy = ZeroFixture(layout), dense = lazy;
	EmptyFields(lazy);
	const double dt = .5 * lazy.StableTimestep();
	Require(lazy.Energy(dt) == 0 && lazy.ElectricEnergy() == 0 && lazy.MagneticEnergy() == 0,
		"initial implicit energy changed");
	EmptyFields(lazy);
	DenseZeros(dense);
	const auto before = dense.HostFieldStorage();
	Require(before.electricBytes >= dense.EdgeCount() * sizeof(double) &&
		before.magneticBytes >= dense.FaceCount() * sizeof(double), "explicit CPU setters failed to materialize");
	for (auto *grid : {&lazy, &dense})
	{
		grid->SetSources({{grid->EdgeCount() / 2, .1}});
	}
	const auto metadata = lazy.DeviceSetupMetadata();
	Require(metadata.electric.empty() && metadata.magnetic.empty(), "metadata retained a global state array");
	EmptyFields(lazy);
	CartesianSubgridDomainPlan plan(lazy, 0, {Base()[0][4]});
	EmptyFields(lazy);
	for (unsigned int step = 0; step < 32; ++step)
	{
		const double pulse = step < 8 ? std::sin(.4 * step) : 0;
		lazy.Step(dt, pulse);
		dense.Step(dt, pulse);
	}
	const auto a = lazy.DeviceStateData(), b = dense.DeviceStateData();
	Require(a.electric == b.electric && a.magnetic == b.magnetic && lazy.Energy(dt) == dense.Energy(dt),
		"CPU lazy fields changed complete-step dynamics");
	std::cout << "IMPLICIT_ZERO_CPU spatial=" << (layout == Grid::CurlLayout::SpatialBlocks)
		<< " exact_fields=1 exact_energy=1 steps=32\n";
}

static void DeviceControl(const SubgridDeviceTarget &target, unsigned int axis, bool domains)
{
	auto lazy = ZeroFixture(Grid::CurlLayout::SpatialBlocks), dense = lazy;
	DenseZeros(dense);
	for (auto *grid : {&lazy, &dense})
	{
		grid->SetSources({{grid->EdgeCount() / 2, .1}});
	}
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend *> pointers;
	for (unsigned int worker = 0; worker < (domains ? 2u : 1u); ++worker)
	{
		target.SelectDomain(worker);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "implicit-state device unavailable");
		target.Verify(*backends.back());
		pointers.push_back(backends.back().get());
	}
	const auto create = [&](const Grid &grid)
	{
		if (domains)
		{
			auto plan = std::make_shared<CartesianSubgridDomainPlan>(grid, axis,
				std::vector<double>{Base()[axis][4]});
			return CreateCartesianSubgridDomains(plan, pointers);
		}
		return backends.front()->CreateCartesianSubgrid(grid);
	};
	auto a = create(lazy), b = create(dense);
	const Grid::Point position{Base()[0][3], Base()[1][3], Base()[2][3]};
	const std::vector<Grid::Integral> functions{lazy.ElectricSample(2, position), lazy.MagneticSample(0, position)};
	CartesianSubgridObservations first(lazy, functions), second(dense, functions);
	first.Bind(*a);
	second.Bind(*b);
	a->ConfigureEnergy(lazy);
	b->ConfigureEnergy(dense);
	EmptyFields(lazy);
	const double dt = .5 * lazy.StableTimestep();
	for (unsigned int step = 0; step < 48; ++step)
	{
		const double pulse = step < 8 ? std::sin(.4 * step) : 0;
		a->Step(dt, pulse);
		b->Step(dt, pulse);
		if (step % 8 == 7)
		{
			std::vector<double> x, y;
			first.Read(x);
			second.Read(y);
			Require(x == y && a->ReadEnergy() == b->ReadEnergy(), "lazy startup changed sparse fields/energy");
			EmptyFields(lazy);
		}
	}
	Require(a->ReadEnergy()[0] + a->ReadEnergy()[1] > 0, "implicit-zero device fixture never excited");
	a->ReadFields(lazy);
	b->ReadFields(dense);
	const auto x = lazy.DeviceStateData(), y = dense.DeviceStateData();
	Require(x.electric == y.electric && x.magnetic == y.magnetic, "explicit full import changed lazy-start fields");
	const auto info = lazy.HostFieldStorage();
	Require(!info.electricImplicitZero && !info.magneticImplicitZero && info.electricBytes > 0 && info.magneticBytes > 0,
		"explicit full import did not materialize host fields");
	std::cout << "IMPLICIT_ZERO_DEVICE axis=" << axis << " domains=" << domains
		<< " sparse_batches=6 exact_fields=1 exact_energy=1 lazy_until_full_import=1\n";
}

int main(int argc, char **argv)
{
	try
	{
		Require(argc <= 2, "unexpected implicit-zero argument");
		const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		StorageControls();
		CPUControl(Grid::CurlLayout::Patterns);
		CPUControl(Grid::CurlLayout::SpatialBlocks);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool domains : {false, true})
			{
				DeviceControl(target, axis, domains);
			}
		}
		std::cout << "PASS: implicit initial fields, exact CPU/device startup and explicit full fallback\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "CSPropMaterial.h"
#include "CSPrimBox.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_mur.h"
#include "subgrid_device_target.h"
#include "tools/constants.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
const std::array<double, 3> epsr{2, 3, 5}, mur{1.2, 1.8, 2.5};

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

template <class Action> static void Reject(Action action, const std::string &expected)
{
	try
	{
		action();
	}
	catch (const std::invalid_argument &error)
	{
		if (std::string(error.what()).find(expected) == std::string::npos)
		{
			throw std::runtime_error("negative control expected " + expected + " but received " + error.what());
		}
		return;
	}
	throw std::runtime_error("unsupported diagonal material was accepted");
}

static Grid::Medium Medium()
{
	Grid::Medium value{epsr[0] * epsilon, mur[0] * mu};
	for (unsigned int n = 0; n < 3; ++n)
	{
		value.epsilonDiagonal[n] = epsr[n] * epsilon;
		value.permeabilityDiagonal[n] = mur[n] * mu;
	}
	return value;
}

static Grid::Axes Lines()
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 6; ++i)
		{
			axes[n].push_back((n + 1) * .001 * i);
		}
	}
	return axes;
}

static std::vector<Grid::Region> RegularRegion()
{
	return {{{{2, 2, 2}, {4, 4, 4}}, {1, 1, 1}}};
}

static void Advance(Grid &grid, double dt, unsigned int steps, const SubgridDeviceTarget *target)
{
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "explicit diagonal material backend missing");
		target->Verify(*backend);
		auto device = backend->CreateCartesianSubgrid(grid);
		Require(bool(device), "explicit backend lacks diagonal subgrid execution");
		for (unsigned int i = 0; i < steps; ++i)
		{
			device->Step(dt, 0);
		}
		device->ReadFields(grid);
	}
	else
	{
		for (unsigned int i = 0; i < steps; ++i)
		{
			grid.Step(dt);
		}
	}
}

//! Exact interior one-step Maxwell curls, independently evaluated in SI units.
static void Yee(const SubgridDeviceTarget *target)
{
	double maximum = 0;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (bool magneticSeed : {false, true})
		{
			const unsigned int a = (n + 1) % 3, b = (n + 2) % 3;
			const auto medium = Medium();
			Grid grid(Lines(), epsilon, mu, RegularRegion(), false,
				[&](const Grid::Point &, const Grid::Point &) { return medium; });
			const auto seed = [&](unsigned int component, const Grid::Point &p)
			{
				return component == b ? p[a] : 0.;
			};
			if (magneticSeed)
			{
				grid.SetMagnetic(seed);
			}
			else
			{
				grid.SetElectric(seed);
			}
			const double dt = .2 * grid.StableTimestep();
			Advance(grid, dt, 1, target);
			const Grid::Point point{.0025, .005, .0075};
			const double actual = magneticSeed ? grid.Electric(n, point) : grid.Magnetic(n, point);
			const double expected = magneticSeed ? dt / medium.Epsilon(n) : -dt / medium.Permeability(n);
			maximum = std::max(maximum, std::abs(actual / expected - 1));
			Require(std::abs(actual / expected - 1) < (target ? 3e-5 : 3e-12),
				"diagonal one-step curl differs from independent Maxwell coefficient");
			Require(std::abs(grid.Evaluate(magneticSeed ? grid.ElectricSample(n, point) :
				grid.MagneticSample(n, point)) / expected - 1) < (target ? 3e-5 : 3e-12),
				"directional point readout and integral weights disagree");
		}
	}
	std::cout << "DIAGONAL_YEE device=" << bool(target) << " directions=6 max_relative=" << maximum << '\n';
}

static void Invariant(const SubgridDeviceTarget *target)
{
	auto medium = [](const Grid::Point &lo, const Grid::Point &)
	{
		auto value = Medium();
		if (lo[0] >= .003)
		{
			value.epsilonDiagonal[1] *= 1.7;
		}
		if (lo[1] >= .006)
		{
			value.permeabilityDiagonal[2] *= 1.4;
		}
		return value;
	};
	const std::vector<Grid::Region> regions{{{{2, 2, 2}, {4, 4, 4}}, {2, 4, 2}}};
	Grid grid(Lines(), epsilon, mu, regions, false, medium);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(371 * p[0] + 537 * p[1] - 891 * p[2] + n);
	});
	const double dt = .6 * grid.StableTimestep(), initial = grid.Energy(dt);
	Advance(grid, dt, 200, target);
	const double drift = std::abs(grid.Energy(dt) / initial - 1);
	Require(drift < (target ? 2e-4 : 2e-11), "diagonal interface modified-energy invariant failed");
	Reject([&] { grid.Step(1.01 * grid.StableTimestep()); }, "timestep");
	std::cout << "DIAGONAL_INTERFACE device=" << bool(target) << " steps=200 drift=" << drift << '\n';
}

static void PaletteAndGuards()
{
	CartesianSubgridMaterialPalette<Grid::Medium> palette;
	palette.Reserve(3);
	auto first = Medium(), second = first, third = first;
	second.epsilonDiagonal[1] = std::nextafter(second.epsilonDiagonal[1], 1.);
	third.permeabilityDiagonal[2] = std::nextafter(third.permeabilityDiagonal[2], 1.);
	for (const auto &value : {first, second, third})
	{
		palette.Append(value);
	}
	palette.Finish();
	Require(palette.UniqueCount() == 3 && palette[1].Epsilon(1) == second.Epsilon(1) &&
		palette[2].Permeability(2) == third.Permeability(2), "directional palette components alias");
	for (unsigned int mode = 0; mode < 4; ++mode)
	{
		auto value = Medium();
		if (mode == 0)
		{
			value.epsilonDiagonal[1] = 0;
		}
		else if (mode == 1)
		{
			value.permeabilityDiagonal[2] = std::numeric_limits<double>::infinity();
		}
		else if (mode == 2)
		{
			value.epsilonDiagonal[0] *= 2;
		}
		else
		{
			value.permeabilityDiagonal[1] = std::numeric_limits<double>::quiet_NaN();
		}
		Reject([&]
		{
			Grid grid(Lines(), epsilon, mu, RegularRegion(), false,
				[&](const Grid::Point &, const Grid::Point &) { return value; });
		}, "diagonal");
	}
	std::cout << "DIAGONAL_PALETTE unique=3 invalid_material_guards=4\n";
}

static void Native(const SubgridDeviceTarget *target, bool magnetic, unsigned int mode = 0)
{
	const std::array<double, 3> electricLoss = mode == 0 ? std::array<double, 3>{.5, .5, .5} :
		mode == 1 ? std::array<double, 3>{0, .3, .9} : std::array<double, 3>{.4, 0, 0};
	const std::array<double, 3> magneticLoss = mode == 0 ? std::array<double, 3>{3000, 3000, 3000} :
		mode == 1 ? std::array<double, 3>{0, 2500, 7000} : std::array<double, 3>{4000, 0, 0};
	// Mode 2 isolates conductivity directionality from epsilon/mu directionality.
	const std::array<double, 3> relativeEpsilon = mode == 2 ? std::array<double, 3>{2, 2, 2} : epsr;
	const std::array<double, 3> relativeMu = mode == 2 ? std::array<double, 3>{1.2, 1.2, 1.2} : mur;
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, i);
		}
	}
	auto *material = new CSPropMaterial(csx.GetParameterSet());
	material->SetIsotropy(false);
	for (unsigned int n = 0; n < 3; ++n)
	{
		material->SetEpsilon(relativeEpsilon[n], n);
		material->SetMue(relativeMu[n], n);
		material->SetKappa(magnetic ? 0 : electricLoss[n], n);
		material->SetSigma(magnetic ? magneticLoss[n] : 0, n);
	}
	csx.AddProperty(material);
	auto *box = new CSPrimBox(csx.GetParameterSet(), material);
	for (unsigned int n = 0; n < 3; ++n)
	{
		box->SetCoord(2 * n, 0.0);
		box->SetCoord(2 * n + 1, 4);
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}));
	op->SetEngineType(target ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "native diagonal geometry rejected");
	int boundary[6];
	std::fill(boundary, boundary + 6, magnetic ? 0 : 1);
	op->SetBoundaryCondition(boundary);
	op->SetTimestepFactor(.4);
	Require(op->CalcECOperator() == 0, "native diagonal operator failed");
	for (std::size_t leaf = 0; leaf < op->Grid().CellCount(); ++leaf)
	{
		const auto value = op->Grid().LeafMedium(leaf);
		for (unsigned int n = 0; n < 3; ++n)
		{
			Require(value.Epsilon(n) == relativeEpsilon[n] * EPS0 &&
				value.Permeability(n) == relativeMu[n] * MUE0,
				"native diagonal CSXCAD/SI material mapping differs");
		}
	}
	const auto seed = [](unsigned int n, const Grid::Point &) { return .2 + .3 * n; };
	if (magnetic)
	{
		op->Grid().SetMagnetic(seed);
	}
	else
	{
		op->Grid().SetElectric(seed);
	}
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "native diagonal explicit backend missing");
		target->Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto *loss = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
	Require(loss, "native diagonal loss provider missing");
	const double dt = op->GetTimestep(), initial = op->Grid().Energy(dt);
	Require(engine->IterateTS(100), "native diagonal evolution failed");
	double error = 0;
	for (unsigned int n = 0; n < 3; ++n)
	{
		const double rate = magnetic ? magneticLoss[n] / (relativeMu[n] * MUE0) :
			electricLoss[n] / (relativeEpsilon[n] * EPS0);
		const double a = .5 * dt * rate, expected = std::pow((1 - a) / (1 + a), 100);
		const Grid::Point point{.0015, .0015, .0015};
		const double value = magnetic ? op->Grid().Magnetic(n, point) : op->Grid().Electric(n, point);
		error = std::max(error, std::abs(value / (.2 + .3 * n) - expected));
	}
	const double balance = std::abs(op->Grid().Energy(dt) + loss->ElectricLoss() +
		loss->MagneticLoss() - initial) / initial;
	Require(error < (target ? 2e-4 : 2e-10) && balance < (target ? 2e-4 : 2e-10),
		"directional constitutive midpoint decay or energy ledger differs");
	std::cout << "DIAGONAL_NATIVE device=" << bool(target) << " magnetic=" << magnetic
		<< " mode=" << mode << " error=" << error << " balance=" << balance << '\n';
}

static void NativeGuards()
{
	for (unsigned int mode = 0; mode < 6; ++mode)
	{
		Reject([&]
		{
			ContinuousStructure csx;
			csx.GetGrid()->SetDeltaUnit(.001);
			for (unsigned int n = 0; n < 3; ++n)
			{
				for (unsigned int i = 0; i <= (n == 0 ? 20u : 4u); ++i)
				{
					csx.GetGrid()->AddDiscLine(n, i);
				}
			}
			auto *material = new CSPropMaterial(csx.GetParameterSet());
			material->SetIsotropy(false);
			for (unsigned int n = 0; n < 3; ++n)
			{
				material->SetEpsilon(n == 1 ? 2 : 1, n);
				material->SetMue(mode == 3 && n == 2 ? -1 : 1, n);
				material->SetKappa(mode == 2 && n == 1 ? -.2 :
					mode == 4 && n == 2 ? std::numeric_limits<double>::infinity() : 0, n);
				material->SetSigma(mode == 5 && n == 1 ? std::numeric_limits<double>::quiet_NaN() : 0, n);
			}
			csx.AddProperty(material);
			auto *box = new CSPrimBox(csx.GetParameterSet(), material);
			for (unsigned int n = 0; n < 3; ++n)
			{
				box->SetCoord(2 * n, 0.0);
				box->SetCoord(2 * n + 1, n == 0 ? 20 : 4);
			}
			Excitation signal;
			signal.SetupDiracPulse(30e9);
			std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({10, 1, 1, 12, 3, 3}));
			op->SetEngineType(Engine::BASIC);
			op->SetExcitationSignal(&signal);
			Require(op->SetGeometryCSX(&csx), "negative native diagonal geometry failed without a reason");
			int boundary[6] = {0, 0, 0, 0, 0, 0};
			if (mode == 0)
			{
				boundary[0] = 3;
				std::array<unsigned int, 6> cells{8, 0, 0, 0, 0, 0};
				op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), cells));
			}
			else if (mode == 1)
			{
				boundary[0] = 2;
				op->AddExtension(new Operator_Ext_CartesianSubgrid_Mur(op.get(), 0));
			}
			op->SetBoundaryCondition(boundary);
			Require(op->CalcECOperator() == 0, "negative native diagonal setup returned an unrelated error");
		}, mode == 0 ? "PML collar" : mode == 1 ? "Mur collar" :
			mode == 3 ? "positive finite diagonal" : "finite nonnegative diagonal");
	}
	std::cout << "DIAGONAL_NATIVE_GUARDS absorber_collars=2 invalid_conductivity=3 invalid_mu=1\n";
}

static void SamplerGuards()
{
	using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
	for (unsigned int mode = 0; mode < 4; ++mode)
	{
		Reject([&]
		{
			ContinuousStructure csx;
			csx.GetGrid()->SetDeltaUnit(.001);
			for (unsigned int n = 0; n < 3; ++n)
			{
				for (unsigned int i = 0; i <= 4; ++i)
				{
					csx.GetGrid()->AddDiscLine(n, i);
				}
			}
			Excitation signal;
			signal.SetupDiracPulse(30e9);
			std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}));
			op->SetExcitationSignal(&signal);
			Require(op->SetGeometryCSX(&csx), "sampler negative-control geometry rejected");
			Loss::Coefficients coefficients{.1, 100};
			coefficients.electricDiagonal = std::array<double, 3>{.1, mode == 0 ? -.2 : .2, .3};
			coefficients.magneticDiagonal = std::array<double, 3>{mode == 2 ? 200. : 100., 200.,
				mode == 1 ? std::numeric_limits<double>::infinity() : 300.};
			if (mode == 3)
			{
				coefficients.electric = std::numeric_limits<double>::quiet_NaN();
			}
			op->AddExtension(new Loss(op.get(), [coefficients](const Grid::Point &, const Grid::Point &)
			{
				return coefficients;
			}));
			Require(op->CalcECOperator() == 0, "sampler negative-control setup failed without a diagnostic");
		}, mode == 2 ? "scalar alias" : "finite nonnegative");
	}
	std::cout << "DIAGONAL_SAMPLER_GUARDS invalid_or_alias=4\n";
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected no argument, --reference, --metal or --cuda");
		SubgridDeviceTarget target;
		if (argc == 2 && std::string(argv[1]) != "--reference")
		{
			target = SubgridDeviceTarget::FromFlag(argv[1]);
		}
		PaletteAndGuards();
		NativeGuards();
		SamplerGuards();
		const std::array<const SubgridDeviceTarget *, 2> devices{nullptr, &target};
		for (const auto *device : devices)
		{
			Yee(device);
			Invariant(device);
			for (unsigned int mode = 0; mode < 3; ++mode)
			{
				Native(device, false, mode);
				Native(device, true, mode);
			}
		}
		std::cout << "PASS: directional Maxwell coefficients, interface energy and native constitutive decay\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}

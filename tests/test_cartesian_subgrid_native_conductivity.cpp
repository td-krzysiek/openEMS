/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "CSPropMaterial.h"
#include "CSPrimBox.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include "tools/constants.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}

//! Uniform E with weak PMC faces: constitutive exponential decay, zero curl.
//! This checks SI kappa mapping and time convergence, not a radiating cavity.
static double Decay(const SubgridDeviceTarget *target, double factor, unsigned int steps,
	double drawingScale = 1, double offset = 0)
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001 / drawingScale);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, (i + offset) * drawingScale);
		}
	}
	constexpr double epsr = 3.5, kappa = .5;
	auto *material = new CSPropMaterial(csx.GetParameterSet());
	material->SetEpsilon(epsr);
	material->SetKappa(kappa);
	csx.AddProperty(material);
	auto *box = new CSPrimBox(csx.GetParameterSet(), material);
	for (unsigned int n = 0; n < 3; ++n)
	{
		box->SetCoord(2 * n, offset * drawingScale);
		box->SetCoord(2 * n + 1, (offset + 4) * drawingScale);
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::vector<double> region(6);
	for (unsigned int n = 0; n < 3; ++n)
	{
		region[n] = (offset + 1) * drawingScale;
		region[n + 3] = (offset + 3) * drawingScale;
	}
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(region));
	op->SetEngineType(target ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "native conductivity geometry rejected");
	int bc[6] = {1, 1, 1, 1, 1, 1};
	op->SetBoundaryCondition(bc);
	op->SetTimestepFactor(factor);
	Require(op->CalcECOperator() == 0, "native conductivity operator failed");
	Require(op->GetNumberOfExtentions() == 1, "native conductivity provider missing/duplicated");
	op->Grid().SetElectric([](unsigned int n, const CartesianSubgrid::Point &) { return .2 + .3 * n; });
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "requested native conductivity device unavailable");
		target->Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto *loss = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
	Require(loss, "native material did not select paired conductivity");
	const double dt = op->GetTimestep(), rate = kappa / (epsr * EPS0);
	const double a = .5 * dt * rate, discrete = std::pow((1 - a) / (1 + a), steps);
	const double continuous = std::exp(-rate * steps * dt), initial = op->Grid().Energy(dt);
	engine->IterateTS(steps / 2);
	engine->IterateTS(steps - steps / 2);
	double error = 0, continuousError = 0;
	for (const auto &point : {CartesianSubgrid::Point{.0005, .0005, .0005},
		CartesianSubgrid::Point{.0015, .0025, .0015}, CartesianSubgrid::Point{.0035, .0035, .0035}})
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double value = op->Grid().Electric(n, point) / (.2 + .3 * n);
			error = std::max(error, std::abs(value - discrete));
			continuousError = std::max(continuousError, std::abs(value - continuous));
		}
	}
	const double balance = std::abs(op->Grid().Energy(dt) + loss->ElectricLoss() + loss->MagneticLoss() - initial) / initial;
	Require(error < (target ? 2e-4 : 2e-10) && balance < (target ? 2e-4 : 2e-10),
		"native conductivity SI mapping or midpoint power differs from analytic decay");
	std::cout << "NATIVE_CONDUCTIVITY device=" << bool(target) << " factor=" << factor << " steps=" << steps
		<< " scale=" << drawingScale << " offset=" << offset << " field_error=" << error
		<< " continuous_error=" << continuousError << " balance=" << balance << std::endl;
	return continuousError;
}


struct MappingFixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	MappingFixture(bool native, bool reverse, double scale, double offset)
	{
		csx.GetGrid()->SetDeltaUnit(.001 / scale);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 4; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, (i + offset) * scale);
			}
		}
		for (unsigned int position = 0; position < 3; ++position)
		{
			const unsigned int index = reverse ? 2 - position : position;
			auto *material = new CSPropMaterial(csx.GetParameterSet());
			material->SetEpsilon(3.5);
			material->SetKappa(native ? index == 0 ? .02 : index == 1 ? .08 : 0 : 0);
			material->SetSigma(native ? index == 0 ? 100 : index == 1 ? 220 : 0 : 0);
			csx.AddProperty(material);
			auto *box = new CSPrimBox(csx.GetParameterSet(), material);
			box->SetPriority(index == 0 ? 1 : index == 1 ? 5 : 9);
			for (unsigned int n = 0; n < 3; ++n)
			{
				const double lo = index == 2 ? 1 : 0;
				const double hi = index == 2 ? (n == 0 ? 1.5 : 2) : index == 1 && n == 0 ? 2 : 4;
				box->SetCoord(2 * n, (lo + offset) * scale);
				box->SetCoord(2 * n + 1, (hi + offset) * scale);
			}
		}
		signal.SetupDiracPulse(30e9);
		std::vector<double> box(6);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box[n] = (offset + 1) * scale;
			box[n + 3] = (offset + 3) * scale;
		}
		op.reset(Operator_CartesianSubgrid::New(box));
		op->SetExcitationSignal(&signal);
		op->SetEngineType(Engine::BASIC);
		Require(op->SetGeometryCSX(&csx), "mapping fixture geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
		if (!native)
		{
			op->AddExtension(new Operator_Ext_CartesianSubgrid_Conductivity(op.get(),
				[](const CartesianSubgrid::Point &lo, const CartesianSubgrid::Point &hi)
				{
					const double x = .5 * (lo[0] + hi[0]), y = .5 * (lo[1] + hi[1]), z = .5 * (lo[2] + hi[2]);
					if (x > .001 && x < .0015 && y > .001 && y < .002 && z > .001 && z < .002)
					{
						return Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{0, 0};
					}
					return x < .002 ? Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{.08, 220} :
						Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{.02, 100};
				}));
		}
		Require(op->CalcECOperator() == 0, "mapping fixture setup failed");
		Require(op->GetNumberOfExtentions() == 1, "native/manual conductivity provider count differs");
	}
	void Seed()
	{
		op->Grid().SetElectric([](unsigned int n, const CartesianSubgrid::Point &p)
		{
			return std::sin(1731 * p[0] + 933 * p[1] - 671 * p[2] + n);
		});
		op->Grid().SetMagnetic([](unsigned int n, const CartesianSubgrid::Point &p)
		{
			return std::cos(633 * p[0] - 731 * p[1] + 411 * p[2] + n) / 376.730313668;
		});
	}
};

static void Mapping(const SubgridDeviceTarget *target, bool reverse, double scale, double offset)
{
	MappingFixture native(true, reverse, scale, offset), expected(false, false, 1, 0);
	// Repeating setup must replace its own provider and retain manual providers.
	Require(native.op->CalcECOperator() == 0 && native.op->GetNumberOfExtentions() == 1,
		"repeated native setup duplicated its conductivity extension");
	native.Seed();
	expected.Seed();
	if (target)
	{
		native.op->SetEngineType(Engine::GPU);
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "requested mapping backend unavailable");
		target->Verify(*backend);
		native.op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> actual(native.op->CreateEngine()), oracle(expected.op->CreateEngine());
	Require(native.op->GetTimestep() == expected.op->GetTimestep(), "drawing units changed the physical timestep");
	const double initial = expected.op->Grid().Energy(expected.op->GetTimestep());
	actual->IterateTS(128);
	oracle->IterateTS(128);
	const auto a = native.op->Grid().DeviceData(), b = expected.op->Grid().DeviceData();
	double error = 0, peak = 0;
	for (const auto pair : {std::make_pair(&a.electric, &b.electric), std::make_pair(&a.magnetic, &b.magnetic)})
	{
		Require(pair.first->size() == pair.second->size(), "drawing units changed mapping dimensions");
		for (std::size_t i = 0; i < pair.first->size(); ++i)
		{
			peak = std::max(peak, std::abs(double(pair.second->at(i))));
			error = std::max(error, std::abs(double(pair.first->at(i)) - pair.second->at(i)));
		}
	}
	auto *loss = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(actual->GetExtension(0));
	auto *expectedLoss = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(oracle->GetExtension(0));
	Require(loss && expectedLoss, "mapping lost its paired extension");
	const double ledger = std::max(std::abs(loss->ElectricLoss() - expectedLoss->ElectricLoss()),
		std::abs(loss->MagneticLoss() - expectedLoss->MagneticLoss())) / initial;
	Require(peak > 0 && error / peak < (target ? 2e-4 : 2e-8) && ledger < (target ? 2e-4 : 2e-8),
		"native kappa/sigma priority or units differ from independent physical sampler");
	std::cout << "NATIVE_MAPPING device=" << bool(target) << " reverse=" << reverse << " scale=" << scale
		<< " offset=" << offset << " field_error=" << error / peak << " ledger_error=" << ledger << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected optional --metal or --cuda");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		const double coarse = Decay(nullptr, .8, 64), fine = Decay(nullptr, .4, 128);
		Require(coarse / fine > 3.8 && coarse / fine < 4.2, "native conductivity time convergence is not second order");
		Decay(nullptr, .8, 64, 1000, 31);
		Decay(&target, .8, 64);
		Decay(&target, .8, 64, 1000, 31);
		Mapping(nullptr, false, 1, 0);
		Mapping(nullptr, true, 1000, 31);
		Mapping(&target, false, 1, 0);
		Mapping(&target, true, 1000, 31);
		std::cout << "PASS: native SI conductivity, analytic decay/time convergence and drawing units" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

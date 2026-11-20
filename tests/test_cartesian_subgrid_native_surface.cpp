/* SPDX-License-Identifier: GPL-3.0-or-later */
#if __has_include("CSPropLossyMetal.h")
#include "subgrid_native_surface_metadata.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropLossyMetal.h"
#include "CSPropMaterial.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;
static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}
struct Result
{
	std::vector<double> electric, magnetic;
	double balance = 0, surfaceLoss = 0, stored = 0;
};
static Result Run(bool bulk, double thickness, const SubgridDeviceTarget *target)
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 2; ++n)
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, i);
		}
	}
	const double copperTop = .8 + thickness / .001;
	for (double z : {0., .5, .8, copperTop, 1.5, 2.})
	{
		csx.GetGrid()->AddDiscLine(2, z);
	}
	auto *metal = new CSPropLossyMetal(csx.GetParameterSet());
	metal->SetConductivity(5.8e7);
	csx.AddProperty(metal);
	auto *box = new CSPrimBox(csx.GetParameterSet(), metal);
	for (unsigned int n = 0; n < 3; ++n)
	{
		box->SetCoord(2 * n, n == 2 ? .8 : 1.5);
		box->SetCoord(2 * n + 1, n == 2 ? copperTop : 2.5);
	}
	box->SetPriority(10);
	if (bulk)
	{
		auto *substrate = new CSPropMaterial(csx.GetParameterSet());
		substrate->SetEpsilon(3.25);
		substrate->SetKappa(.02);
		csx.AddProperty(substrate);
		auto *slab = new CSPrimBox(csx.GetParameterSet(), substrate);
		for (unsigned int n = 0; n < 3; ++n)
		{
			slab->SetCoord(2 * n, n == 2 ? .5 : 0.);
			slab->SetCoord(2 * n + 1, n == 2 ? .8 : 4.);
		}
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({1, 1, .5, 3, 3, 1.5}, {2, 2, 1}));
	op->SetEngineType(target ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "native surface fixture geometry rejected");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	Require(op->CalcECOperator() == 0, "native surface preparation failed");
	// Native Reset consumes the waveform band; reconfigure it for a new setup.
	signal.SetupDiracPulse(30e9);
	Require(op->CalcECOperator() == 0 && op->GetNumberOfExtentions() == 1,
		"repeated native surface setup duplicated its advance provider");
	auto &g = op->Grid();
	Require(g.OpaqueBoundaryFaces().size() == 16, "native surface mapping missed physical body faces");
	std::size_t excluded = 0;
	for (std::size_t leaf = 0; leaf < g.CellCount(); ++leaf)
	{
		excluded += g.LeafIsOpaque(leaf);
	}
	Require(excluded == 4, "native surface mapping changed physical body volume");
	g.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return n == 1 ? .002 * std::exp(-std::pow((p[2] - .0011) / .0004, 2)) : 0;
	});
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "native surface device unavailable");
		target->Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	Engine_Ext_CartesianSubgrid_Conductivity *surface = nullptr;
	unsigned int providers = 0;
	for (std::size_t i = 0; i < engine->GetExtensionCount(); ++i)
	{
		if (auto *provider = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(i)))
		{
			surface = provider;
			++providers;
		}
	}
	Require(providers == 1 && surface, "native surface/bulk forms do not share one advance provider");
	const double dt = op->GetTimestep(), initial = g.Energy(dt);
	Require(initial > 0, "native surface fixture has no incident energy");
	Result result;
	for (unsigned int step = 0; step < 512; ++step)
	{
		engine->IterateTS(1);
		result.balance = std::max(result.balance, std::abs(g.Energy(dt) + surface->AdditionalStoredEnergy() +
			surface->ElectricLoss() + surface->MagneticLoss() + surface->SurfaceLoss() - initial) / initial);
		result.surfaceLoss = surface->SurfaceLoss() / initial;
		result.stored = surface->AdditionalStoredEnergy() / initial;
		if (step % 8 == 0)
		{
			for (std::size_t leaf = 0; leaf < g.CellCount(); ++leaf)
			{
				const auto lo = g.LeafStart(leaf), length = g.LeafLengths(leaf);
				auto centre = lo;
				for (unsigned int n = 0; n < 3; ++n)
				{
					centre[n] += .5 * length[n];
				}
				for (unsigned int n = 0; n < 3; ++n)
				{
					const double e = g.Electric(n, centre), h = g.Magnetic(n, centre);
					if (g.LeafIsOpaque(leaf))
					{
						Require(e == 0 && h == 0, "native metal develops bulk fields");
					}
					result.electric.push_back(e);
					result.magnetic.push_back(h);
				}
			}
		}
	}
	Require(result.balance < (target ? 2e-4 : 2e-9), "native surface/bulk energy ledger failed");
	Require(result.surfaceLoss > 1e-8, "native surface has no measurable Ohmic effect");
	std::cout << "NATIVE_SURFACE bulk=" << bulk << " thickness=" << thickness << " device=" << bool(target) << " excluded=" << excluded
		<< " faces=" << g.OpaqueBoundaryFaces().size() << " balance=" << result.balance
		<< " surface_loss=" << result.surfaceLoss << " stored=" << result.stored << std::endl;
	return result;
}
int main(int argc, char **argv)
{
	if (!SubgridHasNativeSurfaceMetadata<CSPropLossyMetal>::value)
	{
		std::cout << "SKIP: CSXCAD lacks the required volume LossyMetal metadata" << std::endl;
		return 77;
	}

	try
	{
		SetThreadLimit(1);
		const auto target = argc == 2 ? std::string(argv[1]) == "--reference" ? SubgridDeviceTarget{} :
			SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		for (const auto settings : {std::make_pair(false, 18e-6), std::make_pair(true, 18e-6),
			std::make_pair(false, 35e-6), std::make_pair(true, 35e-6)})
		{
			const bool bulk = settings.first;
			const double thickness = settings.second;
			const auto cpu = Run(bulk, thickness, nullptr), device = Run(bulk, thickness, &target);
			double error = 0;
			for (const auto &pair : {std::make_pair(cpu.electric, device.electric), std::make_pair(cpu.magnetic, device.magnetic)})
			{
				Require(pair.first.size() == pair.second.size(), "native surface observation dimensions differ");
				double peak = 0, delta = 0;
				for (std::size_t i = 0; i < pair.first.size(); ++i)
				{
					peak = std::max(peak, std::abs(pair.first[i]));
					delta = std::max(delta, std::abs(pair.first[i] - pair.second[i]));
				}
				error = std::max(error, delta / std::max(peak, 1e-300));
			}
			Require(error < 2e-4, "native surface CPU/device fields disagree");
			std::cout << "NATIVE_SURFACE_COMPARE bulk=" << bulk << " field_error=" << error << std::endl;
		}
		std::cout << "PASS: native smooth volume metal, surface states and shared bulk loss" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

#else
#include <iostream>
int main()
{
	std::cout << "SKIP: CSXCAD has no volume LossyMetal property" << std::endl;
	return 77;
}
#endif

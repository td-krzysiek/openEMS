/* SPDX-License-Identifier: GPL-3.0-or-later */
#if __has_include("CSPropLossyMetal.h")
#include "subgrid_native_surface_metadata.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropLossyMetal.h"
#include "CSPropMaterial.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
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
enum class Shape { Whole, TouchingX, OverlapX, ReverseOverlapX, TouchingZ, PriorityCut, ExplicitCut, DifferentSigma };
struct Result
{
	std::vector<double> fields;
	double loss = 0, balance = 0;
};
static Result Run(Shape shape, const SubgridDeviceTarget *target)
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
	const double lower = .8, upper = .78978, middle = .5 * (lower + upper);
	for (double z : {0., .5, lower, middle, upper, 1.5, 2.})
	{
		csx.GetGrid()->AddDiscLine(2, z);
	}
	auto *metal = new CSPropLossyMetal(csx.GetParameterSet());
	metal->SetConductivity(5.8e7);
	csx.AddProperty(metal);
	const auto addBox = [&](CSProperties *property, double x0, double x1, double z0, double z1, int priority)
	{
		auto *box = new CSPrimBox(csx.GetParameterSet(), property);
		box->SetPriority(priority);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box->SetCoord(2 * n, n == 0 ? x0 : n == 1 ? 1.5 : z0);
			box->SetCoord(2 * n + 1, n == 0 ? x1 : n == 1 ? 2.5 : z1);
		}
	};
	switch (shape)
	{
	case Shape::TouchingX:
	case Shape::DifferentSigma:
		addBox(metal, 1.5, 2, lower, upper, 10);
		if (shape == Shape::DifferentSigma)
		{
			auto *other = new CSPropLossyMetal(csx.GetParameterSet());
			other->SetConductivity(2.9e7);
			csx.AddProperty(other);
			addBox(other, 2, 2.5, lower, upper, 10);
		}
		else
		{
			addBox(metal, 2, 2.5, lower, upper, 10);
		}
		break;
	case Shape::OverlapX:
		addBox(metal, 1.5, 2, lower, upper, 10);
		addBox(metal, 1.5, 2.5, lower, upper, 10);
		break;
	case Shape::ReverseOverlapX:
		addBox(metal, 1.5, 2.5, lower, upper, 10);
		addBox(metal, 1.5, 2, lower, upper, 10);
		break;
	case Shape::TouchingZ:
		addBox(metal, 1.5, 2.5, lower, middle, 10);
		addBox(metal, 1.5, 2.5, middle, upper, 10);
		break;
	case Shape::Whole:
	case Shape::PriorityCut:
	case Shape::ExplicitCut:
		addBox(metal, 1.5, 2.5, shape == Shape::ExplicitCut ? middle : lower, upper, 10);
		break;
	}
	if (shape == Shape::PriorityCut || shape == Shape::ExplicitCut)
	{
		auto *air = new CSPropMaterial(csx.GetParameterSet());
		csx.AddProperty(air);
		addBox(air, 1.5, 2.5, lower, middle, 20);
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({1, 1, .5, 3, 3, 1.5}, {2, 2, 1}));
	op->SetEngineType(target ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "conductor-union native geometry rejected");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	Require(op->CalcECOperator() == 0, "conductor-union preparation failed");
	auto &grid = op->Grid();
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return n == 1 ? .002 * std::exp(-std::pow((p[2] - .0011) / .0004, 2)) : 0;
	});
	if (target)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "conductor-union backend unavailable");
		target->Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto *surface = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
	Require(surface && engine->GetExtensionCount() == 1, "conductor union lost its shared advance provider");
	const double dt = op->GetTimestep(), initial = grid.Energy(dt);
	Require(initial > 0, "conductor-union fixture has no energy");
	Result result;
	for (unsigned int step = 0; step < 1024; ++step)
	{
		engine->IterateTS(1);
		result.loss = (surface->ElectricLoss() + surface->MagneticLoss() + surface->SurfaceLoss()) / initial;
		result.balance = std::max(result.balance,
			std::abs(grid.Energy(dt) + surface->AdditionalStoredEnergy() + surface->ElectricLoss() +
				surface->MagneticLoss() + surface->SurfaceLoss() - initial) / initial);
		if (step % 16 == 0)
		{
			for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
			{
				const auto start = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
				auto centre = start;
				for (unsigned int n = 0; n < 3; ++n)
				{
					centre[n] += .5 * length[n];
				}
				for (unsigned int n = 0; n < 3; ++n)
				{
					result.fields.push_back(grid.Electric(n, centre));
					result.fields.push_back(grid.Magnetic(n, centre) * 376.730313668);
				}
			}
		}
	}
	std::cout << "CONDUCTOR_UNION_STATE shape=" << int(shape) << " device=" << bool(target)
		<< " loss=" << result.loss << " balance=" << result.balance << std::endl;
	Require(result.loss > 1e-8 && result.balance < (target ? 2e-4 : 2e-9),
		"conductor-union loss effect or energy balance failed");
	return result;
}
static void Compare(const Result &a, const Result &b, Shape shape, bool device)
{
	Require(a.fields.size() == b.fields.size(), "equivalent conductor unions changed local dimensions");
	double delta = 0, peak = 0;
	for (std::size_t i = 0; i < a.fields.size(); ++i)
	{
		delta = std::max(delta, std::abs(a.fields[i] - b.fields[i]));
		peak = std::max(peak, std::abs(a.fields[i]));
	}
	const double field = delta / std::max(peak, 1e-300), loss = std::abs(a.loss - b.loss);
	Require(field < (device ? 2e-4 : 2e-10) && loss < (device ? 2e-4 : 2e-10),
		"equivalent native conductor decomposition or priority changed fields/loss");
	std::cout << "CONDUCTOR_UNION shape=" << int(shape) << " device=" << device
		<< " field_error=" << field << " loss_error=" << loss
		<< " balance=" << b.balance << std::endl;
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
		Require(argc <= 2, "expected optional --metal or --cuda");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		const auto whole = Run(Shape::Whole, nullptr);
		for (Shape shape : {Shape::TouchingX, Shape::OverlapX, Shape::ReverseOverlapX, Shape::TouchingZ})
		{
			Compare(whole, Run(shape, nullptr), shape, false);
		}
		const auto cut = Run(Shape::ExplicitCut, nullptr);
		Compare(cut, Run(Shape::PriorityCut, nullptr), Shape::PriorityCut, false);
		Compare(whole, Run(Shape::TouchingZ, &target), Shape::TouchingZ, true);
		Compare(cut, Run(Shape::PriorityCut, &target), Shape::PriorityCut, true);
		bool rejected = false;
		try
		{
			Run(Shape::DifferentSigma, nullptr);
		}
		catch (const std::invalid_argument &error)
		{
			rejected = std::string(error.what()).find("heterogeneous surface conductors") != std::string::npos;
		}
		Require(rejected, "touching heterogeneous conductor stack was silently fitted as one material");
		std::cout << "PASS: native same-metal unions, priority-clipped thickness and heterogeneous-contact guard" << std::endl;
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

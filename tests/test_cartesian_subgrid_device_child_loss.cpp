/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#if __has_include("CSPropLossyMetal.h")
#include "CSPrimBox.h"
#include "CSPropMaterial.h"
#include "CSPropLossyMetal.h"
#include "subgrid_native_surface_metadata.h"
#endif
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;
using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
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
	throw std::runtime_error("invalid constitutive domain operation accepted");
}

static void RequireUnequalChildren(const Grid &grid, unsigned int axes)
{
	bool unequal = false;
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const auto lengths = grid.LeafLengths(leaf);
		for (unsigned int n = 0; n < axes; ++n)
		{
			unequal = unequal || (lengths[n] < .00099 && std::abs(lengths[n] - .0005) > 1e-10);
		}
	}
	Require(unequal, "native loss fixture has no actual unequal local child lengths");
}

struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Engine_Ext_CartesianSubgrid_Conductivity *advance = nullptr;
	Fixture(unsigned int mode, bool ports)
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 4; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, i);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}, Grid::Index{2, 2, 2}));
		Grid::Axes children;
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i < 4; ++i)
			{
				children[n].push_back(i);
				children[n].push_back(i + (i == 1 ? .5 : .43 + .03 * n));
			}
			children[n].push_back(4);
		}
		op->SetChildLines(children);
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "domain loss fixture geometry failed");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
		auto *loss = new Loss(op.get(), [=](const Grid::Point &, const Grid::Point &stop)
		{
			return Loss::Coefficients{mode == 1 || mode == 3 ? (stop[0] <= .002 ? .02 : .06) : 0,
				mode == 2 || mode == 3 ? (stop[1] <= .002 ? 100. : 250.) : 0};
		});
		if (mode >= 3)
		{
			std::vector<Loss::SurfaceAdmittance> surfaces;
			const std::vector<std::pair<Grid::Point, Grid::Point>> paths = {
				{{.001, .0015, .001}, {.0015, .0015, .001}},
				{{.003, .001, .002}, {.004, .001, .002}}};
			for (const auto &path : paths)
			{
				Loss::SurfaceAdmittance surface;
				surface.start = path.first;
				surface.stop = path.second;
				surface.conductance = .01;
				surface.poles = {{10, 1e-10}, {40, 2e-10}, {0, 3e-10}};
				surfaces.push_back(surface);
			}
			loss->SetSurfaceAdmittances(std::move(surfaces));
		}
		op->AddExtension(loss);
		Require(op->CalcECOperator() == 0, "domain loss setup failed");
		auto &g = op->Grid();
		RequireUnequalChildren(g, 3);
		if (ports)
		{
			g.SetResistors({{{.001, .001, .001}, {.001, .001, .002}, 50, 1e-5},
				{{.003, .003, .001}, {.003, .003, .002}, 75, 0}}, op->GetTimestep());
		}
		g.SetElectric([](unsigned int n, const Grid::Point &p)
		{
			return std::sin(731 * p[0] + 423 * p[1] - 511 * p[2] + n);
		});
		g.SetMagnetic([](unsigned int n, const Grid::Point &p)
		{
			return std::cos(933 * p[0] - 331 * p[1] + 211 * p[2] + n) / 376.730313668;
		});
		const auto path = g.VoltagePath({.001, .0015, .001}, {.0015, .0015, .001});
		Require(path.weights.size() > 1, "loss fixture lacks constrained trace");
		// The retained lossless/no-port core does not accumulate soft-source
		// work. Its zero-loss conservation control must therefore be undriven.
		if (mode || ports)
		{
			g.SetSources({{path.weights.front().first, .5}, {path.weights.back().first, -.2},
				{path.weights.front().first, -.1}});
		}
		Require(loss->BuildExtension(), "domain loss port refactor failed");
		engine.reset(op->CreateEngine());
		advance = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
		Require(advance, "domain loss engine pair missing");
	}
};

static double FieldError(const Grid &a, const Grid &b)
{
	const auto x = a.DeviceStateData(), y = b.DeviceStateData();
	double error = 0;
	for (const auto &pair : {std::make_pair(&x.electric, &y.electric), std::make_pair(&x.magnetic, &y.magnetic)})
	{
		double scale = 0, delta = 0;
		Require(pair.first->size() == pair.second->size(), "domain loss field dimensions differ");
		for (std::size_t k = 0; k < pair.first->size(); ++k)
		{
			scale = std::max(scale, std::abs(double(pair.first->at(k))));
			delta = std::max(delta, std::abs(double(pair.first->at(k)) - pair.second->at(k)));
		}
		error = std::max(error, delta / std::max(scale, 1e-300));
	}
	return error;
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis,
	unsigned int mode, bool ports, unsigned int domains, bool mixed)
{
	Fixture cpu(mode, ports), one(mode, ports), split(mode, ports);
	const double dt = cpu.op->GetTimestep();
	const double initial = cpu.op->Grid().Energy(dt);
	const std::vector<double> cuts = domains == 2 ? std::vector<double>{.0021} : std::vector<double>{.0012, .0028};
	if (ports)
	{
		const CartesianSubgridDomainPlan preview(*cpu.advance, axis, cuts);
		const auto owner = std::find_if(preview.Domains().begin(), preview.Domains().end(),
			[](const CartesianSubgridDomainData &d) { return d.fields.rootConductance.empty(); });
		Require(owner != preview.Domains().end(), "no source-only worker for ledger control");
		const auto extra = owner->globalElectric.front();
		for (auto *fixture : {&cpu, &one, &split})
		{
			auto &grid = fixture->op->Grid();
			const auto path = grid.VoltagePath({.001, .0015, .001}, {.0015, .0015, .001});
			grid.SetSources({{path.weights.front().first, .5}, {path.weights.back().first, -.2},
				{path.weights.front().first, -.1}, {extra, .3}});
		}
	}
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(*split.advance, axis, cuts);
	Require((mode ? dt : 0) == plan->RecurrenceTimestep(), "domain loss fixed-step contract missing");
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend*> pointers;
	std::size_t surfaceOwners = 0, sourceOnly = 0;
	for (std::size_t k = 0; k < domains; ++k)
	{
		const auto &d = plan->Domains()[k];
		Require(d.loss.has_value() == (mode != 0), "domain lost constitutive preparation");
		surfaceOwners += d.loss && !d.loss->surfaceRate.empty();
		sourceOnly += d.fields.rootConductance.empty() && !d.fields.sourceEdge.empty();
		target.SelectDomain(k);
		const auto local = mixed && k % 2 ? SubgridDeviceTarget{} : target;
		backends.emplace_back(GPU_Backend::New(local.backend));
		Require(bool(backends.back()), "domain loss backend missing");
		local.Verify(*backends.back());
		pointers.push_back(backends.back().get());
	}
	Require(!ports || sourceOnly > 0, "soft-source work did not reach a separate worker");
	if (mode >= 3 && axis == 0 && !ports)
	{
		Require(surfaceOwners >= 2, "surface histories were not actually distributed");
	}
	auto distributed = CreateCartesianSubgridDomains(plan, pointers);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "single loss backend missing");
	target.Verify(*backend);
	auto single = one.advance->CreateDevice(*backend);
	distributed->RequireFixedTimestep(dt);
	Reject([&] { distributed->Step(dt * .5, 0); });
	Reject([&] { one.advance->PrepareDomainData(); });
	for (unsigned int step = 0; step < 768; ++step)
	{
		const double x = (double(step) - 40) / 12;
		const double drive = std::exp(-x * x);
		cpu.advance->UpdateElectric(dt, drive);
		cpu.advance->UpdateMagnetic(dt);
		single->Step(dt, drive);
		distributed->UpdateElectric(dt, drive);
		if (step == 0)
		{
			Reject([&] { distributed->ReadAdvanceLedger(); });
			Reject([&] { distributed->ReadAdvanceSurfaceLedger(); });
			Reject([&] { distributed->ReadFields(split.op->Grid()); });
		}
		distributed->UpdateMagnetic(dt);
	}
	single->ReadFields(one.op->Grid());
	distributed->ReadFields(split.op->Grid());
	one.advance->ReadDeviceDiagnostics(*single);
	split.advance->ReadDeviceDiagnostics(*distributed);
	Reject([&] { cpu.advance->PrepareDomainData(); });
	Reject([&] { split.advance->PrepareDomainData(); });
	const double singleError = FieldError(one.op->Grid(), split.op->Grid());
	const double cpuError = FieldError(cpu.op->Grid(), split.op->Grid());
	const double scale = std::max(initial, std::abs(cpu.op->Grid().SourceWork()));
	double diagnostics = 0;
	for (const auto &pair : {std::make_pair(cpu.advance, one.advance), std::make_pair(cpu.advance, split.advance)})
	{
		for (const auto &value : {std::make_pair(pair.first->ElectricLoss(), pair.second->ElectricLoss()),
			std::make_pair(pair.first->MagneticLoss(), pair.second->MagneticLoss()),
			std::make_pair(pair.first->SurfaceLoss(), pair.second->SurfaceLoss()),
			std::make_pair(pair.first->AdditionalStoredEnergy(), pair.second->AdditionalStoredEnergy())})
		{
			diagnostics = std::max(diagnostics, std::abs(value.first - value.second) / scale);
		}
	}
	const auto &g = split.op->Grid();
	const double balance = std::abs(g.Energy(dt) + g.ResistorLoss() + split.advance->ElectricLoss() +
		split.advance->MagneticLoss() + split.advance->SurfaceLoss() + split.advance->AdditionalStoredEnergy() -
		g.SourceWork() - initial) / std::max(initial, std::abs(g.SourceWork()));
	std::cout << "CHILD_LOSS_DOMAINS " << target.backend << " axis=" << axis << " mode=" << mode
		<< " ports=" << ports << " domains=" << domains << " mixed=" << mixed << " single_error=" << singleError
		<< " cpu_error=" << cpuError << " ledger_error=" << diagnostics << " balance=" << balance << std::endl;
	Require(singleError < 2e-6 && cpuError < 2e-4 && diagnostics < 2e-4 && balance < 2e-4,
		"constitutive domains field/ledger/passivity gate failed");
	if (mode)
	{
		auto broken = plan->Domains().front();
		broken.loss->magneticRate.push_back(0);
		Reject([&] { CreateCartesianSubgridDomainExecution(*backend, broken); });
		broken = plan->Domains().front();
		if (!broken.loss->conductance.edge.empty())
		{
			broken.loss->conductance.edge.front() = static_cast<std::uint32_t>(broken.ownedElectric);
			Reject([&] { CreateCartesianSubgridDomainExecution(*backend, broken); });
		}
	}
}

#if __has_include("CSPropLossyMetal.h")
struct NativeFixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Engine_Ext_CartesianSubgrid_Conductivity *advance = nullptr;
	NativeFixture(double thickness, bool bulk)
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 2; ++n)
		{
			for (unsigned int i = 0; i <= 4; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, i);
			}
		}
		const double top = .77178 + thickness / .001;
		for (double z : {0., .5, .77178, top, 1.5, 2.})
		{
			csx.GetGrid()->AddDiscLine(2, z);
		}
		auto *metal = new CSPropLossyMetal(csx.GetParameterSet());
		metal->SetConductivity(5.8e7);
		csx.AddProperty(metal);
		auto *box = new CSPrimBox(csx.GetParameterSet(), metal);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box->SetCoord(2 * n, n == 2 ? .77178 : 1.43);
			box->SetCoord(2 * n + 1, n == 2 ? top : 2.57);
		}
		box->SetPriority(10);
		if (bulk)
		{
			auto *substrate = new CSPropMaterial(csx.GetParameterSet());
			substrate->SetEpsilon(3.556);
			substrate->SetKappa(.02);
			csx.AddProperty(substrate);
			auto *slab = new CSPrimBox(csx.GetParameterSet(), substrate);
			for (unsigned int n = 0; n < 3; ++n)
			{
				slab->SetCoord(2 * n, n == 2 ? .5 : 0.);
				slab->SetCoord(2 * n + 1, n == 2 ? .77178 : 4.);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({1, 1, .5, 3, 3, 1.5}, {2, 2, 1}));
		Grid::Axes children;
		for (unsigned int n = 0; n < 2; ++n)
		{
			for (unsigned int i = 0; i < 4; ++i)
			{
				children[n].push_back(i);
				children[n].push_back(i + (i == 2 ? .57 : .43));
			}
			children[n].push_back(4);
		}
		op->SetChildLines(children);
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "native domain metal geometry failed");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
		Require(op->CalcECOperator() == 0 && op->GetNumberOfExtentions() == 1,
			"native domain conductor must prepare exactly one advance");
		RequireUnequalChildren(op->Grid(), 2);
		Require(op->Grid().OpaqueBoundaryFaces().size() == 16, "native domain conductor faces changed");
		op->Grid().SetMagnetic([](unsigned int n, const Grid::Point &p)
		{
			return n == 1 ? .002 * std::exp(-std::pow((p[2] - .0011) / .0004, 2)) : 0;
		});
		engine.reset(op->CreateEngine());
		advance = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
		Require(advance, "native domain conductor advance missing");
	}
};

static void NativeControl(const SubgridDeviceTarget &target, double thickness, bool bulk)
{
	if (!SubgridHasNativeSurfaceMetadata<CSPropLossyMetal>::value)
	{
		std::cout << "SKIP: native domain metal requires optional CSXCAD volume metadata" << std::endl;
		return;
	}
	NativeFixture cpu(thickness, bulk), one(thickness, bulk), split(thickness, bulk);
	const double dt = cpu.op->GetTimestep(), initial = cpu.op->Grid().Energy(dt);
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(*split.advance, 0, std::vector<double>{.002});
	std::size_t poles = 0;
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend*> pointers;
	for (std::size_t k = 0; k < plan->Domains().size(); ++k)
	{
		Require(plan->Domains()[k].loss.has_value(), "native domain conductor dropped its histories");
		poles += plan->Domains()[k].loss->surfaceRate.size();
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "native domain metal backend missing");
		target.Verify(*backends.back());
		pointers.push_back(backends.back().get());
	}
	Require(poles == split.advance->PrepareDomainData().loss->surfaceRate.size() && poles > 0,
		"native constitutive states were duplicated or dropped");
	auto distributed = CreateCartesianSubgridDomains(plan, pointers);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "native single conductor backend missing");
	target.Verify(*backend);
	auto single = one.advance->CreateDevice(*backend);
	for (unsigned int step = 0; step < 512; ++step)
	{
		cpu.advance->UpdateElectric(dt, 0);
		cpu.advance->UpdateMagnetic(dt);
		single->Step(dt, 0);
		distributed->Step(dt, 0);
	}
	single->ReadFields(one.op->Grid());
	distributed->ReadFields(split.op->Grid());
	one.advance->ReadDeviceDiagnostics(*single);
	split.advance->ReadDeviceDiagnostics(*distributed);
	const auto &g = split.op->Grid();
	const double balance = std::abs(g.Energy(dt) + split.advance->ElectricLoss() + split.advance->MagneticLoss() +
		split.advance->SurfaceLoss() + split.advance->AdditionalStoredEnergy() - initial) / initial;
	const double singleError = FieldError(one.op->Grid(), g), cpuError = FieldError(cpu.op->Grid(), g);
	const double lossError = std::abs(cpu.advance->SurfaceLoss() - split.advance->SurfaceLoss()) / initial;
	const double storedError = std::abs(cpu.advance->AdditionalStoredEnergy() - split.advance->AdditionalStoredEnergy()) / initial;
	std::cout << "CHILD_NATIVE_LOSS_DOMAINS " << target.backend << " thickness=" << thickness << " bulk=" << bulk
		<< " poles=" << poles << " single_error=" << singleError << " cpu_error=" << cpuError
		<< " loss_error=" << lossError << " stored_error=" << storedError << " balance=" << balance << std::endl;
	Require(singleError < 2e-6 && cpuError < 2e-4 && lossError < 2e-4 && storedError < 2e-4 && balance < 2e-4,
		"native conductor distributed fields/passivity gate failed");
	Require(split.advance->SurfaceLoss() > initial * 1e-8, "native distributed conductor has no Ohmic effect");
}
#endif

int main(int argc, char **argv)
{
	try
	{
		Require(argc <= 2, "expected one backend flag");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		for (unsigned int mode = 0; mode < 5; ++mode)
		{
			Control(target, mode % 3, mode, false, mode % 2 + 2, false);
			Control(target, (mode + 1) % 3, mode, true, 3, false);
		}
		Control(target, 0, 3, false, 3, false);
		if (target.backend != "reference" && target.physicalDevices == 1)
		{
			Control(target, 0, 3, true, 3, true);
		}
#if __has_include("CSPropLossyMetal.h")
		NativeControl(target, 18e-6, false);
		NativeControl(target, 18e-6, true);
		NativeControl(target, 35e-6, true);
#else
		std::cout << "SKIP: CSXCAD has no volume LossyMetal API" << std::endl;
#endif
		Fixture unsupported(1, false);
		unsupported.op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(unsupported.op.get(), {1, 0, 0, 0, 0, 0}));
		Reject([&] { unsupported.advance->PrepareDomainData(); });
		std::cout << "PASS: native unequal child constitutive states, fields, ports, work and energy" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

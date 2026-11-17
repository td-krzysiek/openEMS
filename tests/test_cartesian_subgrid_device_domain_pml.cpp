/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_energy_test.h"
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
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
	throw std::runtime_error("invalid domain PML operation accepted");
}

struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Engine_Ext_CartesianSubgrid_PML *pml = nullptr;
	Engine_Ext_CartesianSubgrid_Conductivity *advance = nullptr;
	Fixture(unsigned int normal, bool corner, bool loss, bool ports)
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 14; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, i);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({6, 6, 6, 8, 8, 8}, Grid::Index{2, 2, 2}));
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "domain PML geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		std::array<unsigned int, 6> cells{};
		for (unsigned int side = 0; side < 6; ++side)
		{
			if (corner || side / 2 == normal)
			{
				bc[side] = 3;
				cells[side] = 4;
			}
		}
		op->SetBoundaryCondition(bc);
		op->SetTimestepFactor(.8);
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), cells));
		auto *material = loss ? new Loss(op.get(), [](const Grid::Point &lo, const Grid::Point &hi)
		{
			bool inside = true;
			for (unsigned int n = 0; n < 3; ++n)
			{
				const double p = .5 * (lo[n] + hi[n]);
				inside = inside && p >= .006 && p < .008;
			}
			return inside ? Loss::Coefficients{.02, 50.} : Loss::Coefficients{};
		}) : nullptr;
		if (material)
		{
			Loss::SurfaceAdmittance surface;
			surface.start = {.006, .0065, .006};
			surface.stop = {.0065, .0065, .006};
			surface.conductance = .01;
			surface.poles = {{25, 1e-10}, {50, 2e-10}};
			material->SetSurfaceAdmittances({surface});
			op->AddExtension(material);
		}
		Require(op->CalcECOperator() == 0, "domain PML operator failed");
		auto &grid = op->Grid();
		if (ports)
		{
			grid.SetResistors({{{.006, .006, .006}, {.006, .006, .007}, 50, 1e-5},
				{{.008, .008, .006}, {.008, .008, .007}, 75, 0}}, op->GetTimestep());
			Require(material && material->BuildExtension(), "domain PML port preparation failed");
		}
		grid.SetElectric([](unsigned int n, const Grid::Point &p)
		{
			const double x = (p[0] - .007) / .002, y = (p[1] - .007) / .002, z = (p[2] - .007) / .002;
			return (n == 1 ? -z : n == 2 ? y : 0) * std::exp(-x*x - y*y - z*z);
		});
		grid.SetMagnetic([](unsigned int, const Grid::Point &) { return 0.; });
		if (loss)
		{
			const auto path = grid.VoltagePath({.007, .0065, .007}, {.0075, .0065, .007});
			Require(!path.weights.empty(), "domain PML source missing");
			grid.SetSources({{path.weights.front().first, .2}});
		}
		engine.reset(op->CreateEngine());
		for (std::size_t k = 0; k < engine->GetExtensionCount(); ++k)
		{
			if (auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(engine->GetExtension(k)))
			{
				pml = extension;
			}
			if (auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(k)))
			{
				advance = extension;
			}
		}
		Require(pml && bool(advance) == loss, "domain PML engine pair missing");
	}
	void Step(double dt, double drive)
	{
		pml->DoPreVoltageUpdates();
		if (advance)
		{
			advance->UpdateElectric(dt, drive);
		}
		else
		{
			op->Grid().UpdateElectric(dt, drive);
		}
		pml->DoPostVoltageUpdates();
		pml->Apply2Voltages();
		pml->DoPreCurrentUpdates();
		if (advance)
		{
			advance->UpdateMagnetic(dt);
		}
		else
		{
			op->Grid().UpdateMagnetic(dt);
		}
		pml->DoPostCurrentUpdates();
		pml->Apply2Current();
	}
};

static double FieldError(const Grid &a, const Grid &b)
{
	const auto x = a.DeviceStateData(), y = b.DeviceStateData();
	double error = 0;
	for (const auto &pair : {std::make_pair(&x.electric, &y.electric), std::make_pair(&x.magnetic, &y.magnetic)})
	{
		double scale = 0, delta = 0;
		Require(pair.first->size() == pair.second->size(), "domain PML field dimensions differ");
		for (std::size_t k = 0; k < pair.first->size(); ++k)
		{
			scale = std::max(scale, std::abs(double(pair.first->at(k))));
			delta = std::max(delta, std::abs(double(pair.first->at(k)) - pair.second->at(k)));
		}
		error = std::max(error, delta / std::max(scale, 1e-300));
	}
	return error;
}

static void Control(const SubgridDeviceTarget &target, unsigned int normal, unsigned int axis,
	bool corner, bool loss, bool ports, unsigned int domains, bool mixed)
{
	Fixture cpu(normal, corner, loss, ports), one(normal, corner, loss, ports), split(normal, corner, loss, ports);
	const double dt = cpu.op->GetTimestep(), initial = cpu.op->Grid().Energy(dt);
	const auto cuts = domains == 2 ? std::vector<double>{.0061} : std::vector<double>{.0045, .0081};
	auto plan = loss ? std::make_shared<CartesianSubgridDomainPlan>(*split.advance, *split.pml, axis, cuts) :
		std::make_shared<CartesianSubgridDomainPlan>(*split.pml, axis, cuts);
	Require(plan->PMLTimestep() == dt, "domain PML timestep missing");
	std::size_t ownedE = 0, ownedH = 0;
	const auto whole = loss ? split.pml->PrepareDomainData(split.advance) : split.pml->PrepareDomainData();
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend*> devices;
	for (std::size_t k = 0; k < domains; ++k)
	{
		const auto &domain = plan->Domains()[k];
		Require(domain.pml.has_value(), "domain PML histories lost");
		ownedE += domain.pml->electric.target.size();
		ownedH += domain.pml->magnetic.target.size();
		target.SelectDomain(k);
		const auto local = mixed && k % 2 ? SubgridDeviceTarget{} : target;
		backends.emplace_back(GPU_Backend::New(local.backend));
		Require(bool(backends.back()), "domain PML backend unavailable");
		local.Verify(*backends.back());
		std::cout << "PML_DOMAIN_BACKEND owner=" << k << " backend=" << backends.back()->GetName() << std::endl;
		devices.push_back(backends.back().get());
	}
	Require(ownedE == whole.electric.target.size() && ownedH == whole.magnetic.target.size(),
		"PML histories were duplicated or lost");
	auto distributed = CreateCartesianSubgridDomains(plan, devices);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "single PML backend missing");
	target.Verify(*backend);
	std::cout << "PML_SINGLE_BACKEND backend=" << backend->GetName() << std::endl;
	auto single = loss ? one.advance->CreateDevice(*backend) : backend->CreateCartesianSubgrid(one.op->Grid());
	single->ConfigureEnergy(one.op->Grid());
	distributed->ConfigureEnergy(split.op->Grid());
	single->RequireFixedTimestep(dt);
	// Deliberately leave another context current before allocation, every
	// hook and diagnostics. The PML component must select its own owner.
	backends.back()->Activate();
	auto pml = one.pml->CreateDevice(*backend, single->Fields());
	const auto hook = [&](bool magnetic, GPU_Hook phase, unsigned int step)
	{
		backends.back()->Activate();
		pml->Hook(magnetic, phase, Engine_Box{}, step);
	};
	distributed->RequireFixedTimestep(dt);
	Reject([&] { distributed->Step(dt * .5, 0); });
	Reject([&] { one.pml->PrepareDomainData(one.advance); });
	if (loss)
	{
		Reject([&] { split.advance->PrepareDomainData(); });
		Reject([&] { split.pml->PrepareDomainData(); });
	}
	for (unsigned int step = 0; step < 512; ++step)
	{
		const double drive = loss ? std::exp(-std::pow((double(step) - 20) / 8, 2)) : 0;
		cpu.Step(dt, drive);
		hook(false, GPU_PRE, step);
		single->UpdateElectric(dt, drive);
		hook(false, GPU_POST, step);
		hook(false, GPU_APPLY, step);
		hook(true, GPU_PRE, step);
		single->UpdateMagnetic(dt);
		hook(true, GPU_POST, step);
		hook(true, GPU_APPLY, step);
		distributed->UpdateElectric(dt, drive);
		if (step == 0)
		{
			Reject([&] { distributed->ReadBoundaryWorkOut(); });
		}
		distributed->UpdateMagnetic(dt);
	}
	single->ReadFields(one.op->Grid());
	distributed->ReadFields(split.op->Grid());
	CheckSubgridEnergy(*single, one.op->Grid(), 1, "pml-single");
	CheckSubgridEnergy(*distributed, split.op->Grid(), static_cast<unsigned int>(devices.size()), "pml-domains");
	const double singleError = FieldError(one.op->Grid(), split.op->Grid());
	const double cpuError = FieldError(cpu.op->Grid(), split.op->Grid());
	const double work = distributed->ReadBoundaryWorkOut();
	const double scale = std::max(initial, std::abs(cpu.op->Grid().SourceWork()));
	backends.back()->Activate();
	const double workError = std::abs(work - pml->ReadBoundaryWorkOut()) / scale;
	const double cpuWorkError = std::abs(work - cpu.pml->BoundaryWorkOut()) / scale;
	double material = 0, stored = 0;
	if (loss)
	{
		one.advance->ReadDeviceDiagnostics(*single);
		split.advance->ReadDeviceDiagnostics(*distributed);
		material = split.advance->ElectricLoss() + split.advance->MagneticLoss() + split.advance->SurfaceLoss();
		stored = split.advance->AdditionalStoredEnergy();
		Require(material > 0 && split.advance->SurfaceLoss() > 0, "distributed material loss has no effect");
	}
	const auto &g = split.op->Grid();
	const double balance = std::abs(g.Energy(dt) + stored + material + work + g.ResistorLoss() -
		g.SourceWork() - initial) / scale;
	Require(singleError < 2e-6 && cpuError < 2e-4 && workError < 2e-4 && cpuWorkError < 2e-4 && balance < 2e-4,
		"domain PML fields, signed work or energy balance changed");
	Require(std::abs(work) > initial * 1e-3, "PML boundary work has no effect");
	auto broken = plan->Domains().front();
	Require(!broken.pml->electric.target.empty(), "negative domain PML fixture lacks target");
	broken.pml->electric.target.front() = static_cast<std::uint32_t>(broken.ownedElectric);
	Reject([&] { CreateCartesianSubgridWorker(*backend, broken); });
	Reject([&] { cpu.pml->PrepareDomainData(cpu.advance); });
	std::cout << "PML_DOMAINS normal=" << normal << " axis=" << axis << " corner=" << corner
		<< " loss=" << loss << " ports=" << ports << " domains=" << domains << " mixed=" << mixed
		<< " electric=" << ownedE << " magnetic=" << ownedH << " single_error=" << singleError
		<< " cpu_error=" << cpuError << " work_error=" << workError << " cpu_work_error=" << cpuWorkError
		<< " balance=" << balance
		<< std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		Require(argc <= 2, "unexpected domain PML arguments");
		if (target.physicalDevices > 1)
		{
			for (unsigned int axis = 0; axis < 3; ++axis)
			{
				Control(target, axis, axis, true, false, false, 2, false);
			}
			Control(target, 0, 0, true, true, true, 2, false);
			Control(target, 1, 2, true, true, true, 2, false);
		}
		else
		{
			for (unsigned int axis = 0; axis < 3; ++axis)
			{
				Control(target, axis, (axis + 1) % 3, false, false, false, 2, false);
				Control(target, axis, axis, true, true, true, 3, false);
			}
		}
		if (target.backend != "reference" && target.physicalDevices == 1)
		{
			Control(target, 0, 0, true, true, true, 3, true);
		}
		Fixture zero(0, false, false, false), other(0, false, true, true);
		Reject([&] { CartesianSubgridDomainPlan invalid(*other.advance, *zero.pml, 0, {.0061}); });
		zero.op->Grid().SetElectric([](unsigned int, const Grid::Point &) { return 0.; });
		zero.op->Grid().SetMagnetic([](unsigned int, const Grid::Point &) { return 0.; });
		zero.Step(zero.op->GetTimestep(), 0);
		Require(zero.pml->BoundaryWorkOut() == 0, "zero-state PML control unexpectedly acquired work");
		Reject([&] { zero.pml->PrepareDomainData(); });
		std::cout << "PASS_REJECT_WARM_ZERO_PML_AND_MISMATCHED_GRID" << std::endl;
		std::cout << "PASS: target-owned PML histories, ordered halos and boundary work" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

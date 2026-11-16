/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid_runtime.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include "openems.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using Grid = CartesianSubgrid;
static void Require(bool value, const char *message)
{
	if (!value)
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
	throw std::runtime_error("invalid native domain option was accepted");
}

static void Options()
{
	using Config = CartesianSubgridRuntime;
	const auto value = Config::Parse("y", "-2,3.5", "auto@1,reference,auto@0");
	Require(value.axis == 1 && value.drawingCuts == std::vector<double>{-2, 3.5} &&
		value.workers[0].device == 1 && value.workers[1].backend == "reference", "domain parsing changed");
	Reject([&] { value.CreateBackends(true); });
	for (const auto &cuts : {"0,0", "1,0", "nan", "inf", "1junk", "1,", " 1"})
	{
		Reject([&] { Config::Parse("x", cuts, "reference,reference"); });
	}
	for (const auto &workers : {"", "cpu,reference", "reference@0,reference", "auto@-1,reference",
		"auto@1junk,reference", "auto@2147483648,reference", "auto@,reference", "auto,,reference"})
	{
		Reject([&] { Config::Parse("x", "1", workers); });
	}
	Reject([] { Config::Parse("X", "1", "reference,reference"); });
	Require(Config::Parse("z", "", "reference").workers.size() == 1, "single explicit domain rejected");
	Require(Config::Parse("z", "none", "reference").drawingCuts.empty(), "explicit no-cut sentinel rejected");
	Reject([] { Config::Parse("x", "none", "reference,reference"); });
	openEMS settings;
	settings.SetLibraryArguments({"subgrid-domain-axis=x"});
	Reject([] { Config::ReadOptions(); });
	settings.SetLibraryArguments({"subgrid-domain-axis=z", "subgrid-domain-cuts=-2,3.5",
		"subgrid-domain-workers=reference,reference,reference"});
	Require(Config::ReadOptions()->axis == 2, "library domain options missing");
	settings.SetLibraryArguments({"subgrid-domain-axis=x", "subgrid-domain-cuts=none",
		"subgrid-domain-workers=reference"});
	Require(Config::ReadOptions()->workers.size() == 1 && Config::ReadOptions()->drawingCuts.empty(),
		"one-worker library sentinel rejected");
	settings.SetLibraryArguments({});
	Require(!Config::ReadOptions(), "domain settings leaked to the next library run");
}

struct Fixture
{
	ContinuousStructure geometry;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Engine_Ext_CartesianSubgrid_PML *pml = nullptr;
	Fixture(const SubgridDeviceTarget &target, unsigned int axis, bool loss,
		double unit, double offset, bool domains)
	{
		geometry.GetGrid()->SetDeltaUnit(unit);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 14; ++i)
			{
				geometry.GetGrid()->AddDiscLine(n, offset + i);
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({offset+6, offset+6, offset+6,
			offset+8, offset+8, offset+8}));
		op->SetEngineType(Engine::GPU);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&geometry), "native domain geometry rejected");
		int boundaries[6] = {3, 3, 3, 3, 3, 3};
		op->SetBoundaryCondition(boundaries);
		op->SetTimestepFactor(.8);
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), {4, 4, 4, 4, 4, 4}));
		if (loss)
		{
			op->AddExtension(new Operator_Ext_CartesianSubgrid_Conductivity(op.get(),
				[unit](const Grid::Point &lo, const Grid::Point &hi)
				{
					bool inside = true;
					for (unsigned int n = 0; n < 3; ++n)
					{
						const double p = .5 * (lo[n] + hi[n]) / unit;
						inside = inside && p >= 6 && p < 8;
					}
					return inside ? Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{.02, 50} :
						Operator_Ext_CartesianSubgrid_Conductivity::Coefficients{};
				}));
		}
		Require(op->CalcECOperator() == 0, "native domain operator failed");
		signal.SetupDiracPulse(30e9);
		Require(signal.buildExcitationSignal(1024), "native domain waveform failed");
		op->Grid().SetElectric([unit](unsigned int n, const Grid::Point &p)
		{
			return (n == 1 ? 1. : -.5) * std::exp(-std::pow((p[0]/unit-7)/2, 2) -
				std::pow((p[1]/unit-7)/2, 2) - std::pow((p[2]/unit-7)/2, 2));
		});
		if (domains)
		{
			std::vector<std::unique_ptr<GPU_Backend>> backends;
			for (std::size_t k = 0; k < 2; ++k)
			{
				target.SelectDomain(k);
				backends.emplace_back(GPU_Backend::New(target.backend));
				Require(bool(backends.back()), "native worker unavailable");
				target.Verify(*backends.back());
			}
			op->SetDomainExecution(axis, {offset+6.1}, std::move(backends));
		}
		else
		{
			target.SelectDomain(0);
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
			Require(bool(backend), "native single backend unavailable");
			target.Verify(*backend);
			op->SetGPUBackend(backend.release());
		}
		engine.reset(op->CreateEngine());
		for (std::size_t k = 0; k < engine->GetExtensionCount(); ++k)
		{
			if (auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(engine->GetExtension(k)))
			{
				pml = extension;
			}
		}
		Require(pml, "native PML mirror missing");
	}
};

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool loss)
{
	const double offset = axis == 1 ? 20. : -7.;
	const double unit = axis == 2 ? 1e-6 : 1e-3;
	Fixture one(target, axis, loss, unit, offset, false), split(target, axis, loss, unit, offset, true);
	for (unsigned int batch = 0; batch < 8; ++batch)
	{
		one.engine->IterateTS(32);
		split.engine->IterateTS(32);
		const auto a = one.op->Grid().DeviceStateData(), b = split.op->ObservationGrid().DeviceStateData();
		Require(a.electric == b.electric && a.magnetic == b.magnetic, "native complete-batch fields differ");
		Require(one.engine->GetNumberOfTimesteps() == split.engine->GetNumberOfTimesteps(), "native epoch differs");
		double e = 0, h = 0, es = 0, hs = 0;
		const unsigned int stride[3] = {1, 1, 1};
		one.engine->CalcFastEnergy(stride, e, h);
		split.engine->CalcFastEnergy(stride, es, hs);
		Require(e == es && h == hs, "native complete-batch energy differs");
		Require(std::abs(one.pml->BoundaryWorkOut() - split.pml->BoundaryWorkOut()) /
			std::max(std::abs(one.pml->BoundaryWorkOut()), 1e-300) < 2e-4,
			"native PML boundary-work mirror differs");
	}
	Require(std::abs(split.pml->BoundaryWorkOut()) > 0, "native PML work was never read back");
	split.op->InvalidateObservations();
	Reject([&] { split.op->ObservationGrid(); });
	std::cout << "NATIVE_DOMAINS axis=" << axis << " loss=" << loss << " unit=" << unit
		<< " offset=" << offset << " steps=256 exact_fields=1 exact_energy=1" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		Options();
		Require(argc <= 2, "unexpected native-runtime arguments");
		const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool loss : {false, true})
			{
				Control(target, axis, loss);
			}
		}
		std::cout << "PASS: native fixed domains, drawing-coordinate cuts, PML/loss mirrors and options" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

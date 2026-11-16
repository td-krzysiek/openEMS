/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "tools/useful.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}
template <class F> static void Reject(F call)
{
	try
	{
		call();
	}
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid conductivity input/phase was admitted");
}
struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	Fixture()
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
		op.reset(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}));
		op->SetExcitationSignal(&signal);
		op->SetEngineType(Engine::BASIC);
		Require(op->SetGeometryCSX(&csx), "guard fixture geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
		Require(op->CalcECOperator() == 0, "guard fixture setup rejected");
	}
};
static Loss::Sampler Constant(double electric, double magnetic)
{
	return [=](const CartesianSubgrid::Point &, const CartesianSubgrid::Point &)
	{
		return Loss::Coefficients{electric, magnetic};
	};
}
int main()
{
	try
	{
		SetThreadLimit(1);
		Fixture f;
		Reject([&] { Loss missing(nullptr, Constant(0, 0)); });
		Reject([&] { Loss missing(f.op.get(), {}); });
		for (double value : {-1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
		{
			Reject([&] { Loss bad(f.op.get(), Constant(value, 0)); bad.BuildExtension(); });
			Reject([&] { Loss bad(f.op.get(), Constant(0, value)); bad.BuildExtension(); });
		}
		Reject([&] { Loss overflow(f.op.get(), Constant(0, 1e308)); overflow.BuildExtension(); });
		f.op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(f.op.get(), {0, 1, 0, 0, 0, 0}));
		Reject([&] { Loss history(f.op.get(), Constant(.02, 0)); history.BuildExtension(); });
		std::cout << "PASS: sampler/sign/finite/rate/absorber-combination guards" << std::endl;

		Fixture valid;
		auto *loss = new Loss(valid.op.get(), Constant(.02, 100));
		valid.op->AddExtension(loss);
		Require(loss->BuildExtension(), "valid conductivity setup failed");
		std::unique_ptr<Engine> engine(valid.op->CreateEngine());
		auto *pair = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(0));
		Require(pair, "paired conductivity engine missing");
		const double dt = valid.op->GetTimestep();
		Reject([&] { pair->UpdateElectric(dt / 2, 0); });
		Reject([&] { pair->UpdateMagnetic(dt); });
		pair->UpdateElectric(dt, 0);
		Reject([&] { pair->UpdateElectric(dt, 0); });
		Reject([&] { pair->UpdateMagnetic(dt / 2); });
		pair->UpdateMagnetic(dt);

		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New("reference"));
		Require(bool(backend), "reference backend unavailable");
		auto device = pair->CreateDevice(*backend);
		Reject([&] { device->UpdateElectric(dt / 2, 0); });
		Reject([&] { device->UpdateMagnetic(dt); });
		device->UpdateElectric(dt, 0);
		Reject([&] { device->ReadFields(valid.op->Grid()); });
		Reject([&] { device->ReadAdvanceLedger(); });
		Reject([&] { device->UpdateElectric(dt, 0); });
		device->UpdateMagnetic(dt);
		device->ReadFields(valid.op->Grid());
		const auto before = device->ReadAdvanceLedger();
		Require(before == device->ReadAdvanceLedger(), "repeated observation changed a device loss ledger");
		std::cout << "PASS: CPU/reference fixed-step, phase and repeated-observation guards" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

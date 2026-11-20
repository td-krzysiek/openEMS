/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

static void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::runtime_error(message);
	}
}

static void Control(Excitation::ExciteTypes kind, bool device, const SubgridDeviceTarget &target)
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
	op->SetEngineType(device ? Engine::GPU : Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "waveform geometry setup failed");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	Require(op->CalcECOperator() == 0, "waveform operator setup failed");
	if (kind == Excitation::Step)
	{
		signal.SetupStepExcite(30e9);
	}
	else if (kind == Excitation::DiracPulse)
	{
		signal.SetupDiracPulse(30e9);
	}
	else if (kind == Excitation::CustomExcite)
	{
		signal.SetupCustomExcite("1+t*1e10", 30e9, 30e9);
	}
	else
	{
		signal.SetupSinusoidal(30e9);
	}
	Require(signal.buildExcitationSignal(4) && signal.GetLength() > 0,
		"waveform control requires a real stored signal");
	if (kind == Excitation::Step || kind == Excitation::DiracPulse)
	{
		Require(signal.GetLength() == 2, "pulse/Step representation changed");
	}
	op->Grid().SetSources({{op->Grid().EdgeCount() / 2, .01}});
	auto expected = op->Grid(), legacy = expected;
	const auto period = signal.GetSignalPeriodTS();
	const unsigned int steps = period ? 4 * period + 3 : 64;
	Require(steps > signal.GetLength(), "waveform test must advance past stored samples");
	for (unsigned int i = 0; i < steps; ++i)
	{
		const auto position = period ? i % period : i;
		double drive = 0;
		if (kind == Excitation::Step)
		{
			drive = 1;
		}
		else if (kind == Excitation::DiracPulse)
		{
			drive = i == 1 ? 1 : 0;
		}
		else
		{
			drive = signal.GetVoltageSignal()[position < signal.GetLength() ? position : 0];
		}
		expected.Step(op->GetTimestep(), drive);
		legacy.Step(op->GetTimestep(), i < signal.GetLength() ? signal.GetVoltageSignal()[i] : 0);
	}
	if (device)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
		Require(bool(backend), "requested waveform backend unavailable");
		target.Verify(*backend);
		op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	engine->IterateTS(3);
	engine->IterateTS(steps - 3);
	const auto a = expected.DeviceData(), b = op->Grid().DeviceData(), old = legacy.DeviceData();
	double scale = 0, error = 0, legacyError = 0;
	for (const auto pair : {std::make_pair(&a.electric, &b.electric), std::make_pair(&a.magnetic, &b.magnetic)})
	{
		for (std::size_t i = 0; i < pair.first->size(); ++i)
		{
			scale = std::max(scale, std::abs(double(pair.first->at(i))));
			error = std::max(error, std::abs(double(pair.first->at(i)) - pair.second->at(i)));
		}
	}
	for (std::size_t i = 0; i < a.electric.size(); ++i)
	{
		legacyError = std::max(legacyError, std::abs(double(a.electric[i]) - old.electric[i]));
	}
	Require(scale > 0 && error / scale < 2e-4, "native waveform differs from intended core drive");
	if (kind != Excitation::DiracPulse)
	{
		Require(legacyError / scale > 1e-3, "waveform control cannot distinguish the legacy zero-tail bug");
	}
	else
	{
		Require(legacyError == 0, "Dirac control should retain the existing zero-tail evolution");
	}
	std::cout << "WAVEFORM type=" << kind << " device=" << device << " period=" << period
		<< " stored=" << signal.GetLength() << " steps=" << steps << " field_error=" << error / scale
		<< " legacy_error=" << legacyError / scale << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected optional --metal or --cuda");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		for (const auto kind : {Excitation::Step, Excitation::DiracPulse, Excitation::CustomExcite, Excitation::Sinusoidal})
		{
			Control(kind, false, target);
			Control(kind, true, target);
		}
		std::cout << "PASS: native sample-zero tails, periodic sources and unchanged Dirac control" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

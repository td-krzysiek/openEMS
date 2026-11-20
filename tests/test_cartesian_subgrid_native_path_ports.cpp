/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropLumpedElement.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}
struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	Fixture(bool native, bool reverse = false, bool partial = false)
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 4; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, i);
			}
		}
		if (native)
		{
			for (unsigned int i = 0; i < 2; ++i)
			{
				auto *resistor = new CSPropLumpedElement(csx.GetParameterSet());
				resistor->SetResistance(i ? 75 : 50);
				resistor->SetDirection(2);
				resistor->SetLEtype(CSPropLumpedElement::PARALLEL);
				csx.AddProperty(resistor);
				auto *box = new CSPrimBox(csx.GetParameterSet(), resistor);
				const double lo = i ? 1 : 0, hi = i ? (partial ? 2.1 : 3) : 4;
				for (unsigned int n = 0; n < 3; ++n)
				{
					box->SetCoord(2 * n, n == 2 ? (reverse ? hi : lo) : n == 0 ? 2 + i : 2);
					box->SetCoord(2 * n + 1, n == 2 ? (reverse ? lo : hi) : n == 0 ? 2 + i : 2);
				}
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New({1, 1, 1, 3, 3, 3}));
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "native path-port geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
		Require(op->CalcECOperator() == 0, "native path-port preparation failed");
		if (!native)
		{
			std::vector<Grid::Resistor> cells;
			const std::vector<double> first{0, .001, .0015, .002, .0025, .003, .004};
			for (std::size_t i = 1; i < first.size(); ++i)
			{
				cells.push_back({{.002, .002, first[i - 1]}, {.002, .002, first[i]}, 50 * (first[i] - first[i - 1]) / .004, 0});
			}
			cells.push_back({{.003, .002, .001}, {.003, .002, .002}, 75. / 2, 0});
			cells.push_back({{.003, .002, .002}, {.003, .002, .003}, 75. / 2, 0});
			op->Grid().SetResistors(cells, op->GetTimestep());
		}
		op->Grid().SetElectric([](unsigned int n, const Grid::Point &p)
		{
			return std::sin(1731 * p[0] + 933 * p[1] - 671 * p[2] + n);
		});
		op->Grid().SetMagnetic([](unsigned int n, const Grid::Point &p)
		{
			return std::cos(633 * p[0] - 731 * p[1] + 411 * p[2] + n) / 376.730313668;
		});
	}
};
static void Compare(bool reverse, const SubgridDeviceTarget *target)
{
	Fixture native(true, reverse), manual(false);
	if (target)
	{
		native.op->SetEngineType(Engine::GPU);
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
		Require(bool(backend), "path-port backend unavailable");
		target->Verify(*backend);
		native.op->SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> actual(native.op->CreateEngine()), expected(manual.op->CreateEngine());
	const double dt = native.op->GetTimestep(), initial = native.op->Grid().Energy(dt);
	Require(dt == manual.op->GetTimestep() && initial > 0, "path-port reference timestep/energy differs");
	double balance = 0, delta = 0, peak = 0;
	for (unsigned int step = 0; step < 1024; ++step)
	{
		actual->IterateTS(1);
		expected->IterateTS(1);
		const auto &g = native.op->Grid();
		balance = std::max(balance, std::abs(g.Energy(dt) + g.ResistorLoss() - initial) / initial);
		for (unsigned int port = 0; port < 2; ++port)
		{
			const Grid::Point a{port ? .003 : .002, .002, port ? .001 : 0};
			const Grid::Point b{port ? .003 : .002, .002, port ? .003 : .004};
			const double actualVoltage = g.Evaluate(g.VoltagePath(reverse ? b : a, reverse ? a : b));
			const double referenceVoltage = manual.op->Grid().Evaluate(manual.op->Grid().VoltagePath(a, b));
			const double wanted = (reverse ? -1 : 1) * referenceVoltage;
			delta = std::max(delta, std::abs(actualVoltage - wanted));
			peak = std::max(peak, std::abs(wanted));
		}
	}
	const double fieldError = delta / std::max(peak, 1e-300);
	const double lossError = std::abs(native.op->Grid().ResistorLoss() - manual.op->Grid().ResistorLoss()) / initial;
	Require(fieldError < (target ? 2e-4 : 2e-10) && lossError < (target ? 2e-4 : 2e-10) &&
		balance < (target ? 2e-4 : 2e-10), "native multi-edge port voltage, loss or energy balance differs");
	std::cout << "NATIVE_PATH_PORT reverse=" << reverse << " device=" << bool(target)
		<< " voltage_error=" << fieldError << " loss_error=" << lossError << " balance=" << balance << std::endl;
}
int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected optional --metal or --cuda");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		for (bool reverse : {false, true})
		{
			Compare(reverse, nullptr);
			Compare(reverse, &target);
		}
		bool rejected = false;
		try
		{
			Fixture partial(true, false, true);
		}
		catch (const std::invalid_argument &error)
		{
			rejected = std::string(error.what()).find("complete local-edge path") != std::string::npos;
		}
		Require(rejected, "partial-edge native resistor was admitted");
		std::cout << "PASS: native multi-edge resistor paths, orientations, manual physical-length series-R oracle and partial-edge guard" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

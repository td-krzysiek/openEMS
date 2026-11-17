/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPropLumpedElement.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
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
	throw std::runtime_error("invalid rectangular port geometry was accepted");
}
static Grid::Axes Lines(bool graded, bool fine)
{
	Grid::Axes result;
	for (unsigned int n = 0; n < 3; ++n)
	{
		const double base[5] = {0, graded ? .7 + .1 * n : 1, graded ? 1.8 + .1 * n : 2,
			graded ? 3.1 + .1 * n : 3, 4};
		for (unsigned int k = 0; k < 4; ++k)
		{
			result[n].push_back(base[k]);
			if (fine)
			{
				result[n].push_back(.5 * (base[k] + base[k + 1]));
			}
		}
		result[n].push_back(4);
	}
	return result;
}
struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator> op;
	std::unique_ptr<Engine> engine;
	Grid::Axes fine;
	Grid::Point start{}, stop{};
	Fixture(bool compact, bool graded, unsigned int direction, bool sheet, bool caps, bool reverse, bool refined)
		: fine(Lines(graded, true))
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		const auto axes = Lines(graded, !compact || !refined);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (double value : axes[n])
			{
				csx.GetGrid()->AddDiscLine(n, value);
			}
			start[n] = fine[n][refined ? 3 : 1];
			stop[n] = fine[n][refined ? 5 : 7];
		}
		if (sheet)
		{
			start[(direction + 2) % 3] = stop[(direction + 2) % 3] = fine[(direction + 2) % 3][4];
		}
		auto *resistor = new CSPropLumpedElement(csx.GetParameterSet());
		resistor->SetResistance(50);
		resistor->SetDirection(direction);
		resistor->SetCaps(caps);
		resistor->SetLEtype(CSPropLumpedElement::PARALLEL);
		csx.AddProperty(resistor);
		auto *box = new CSPrimBox(csx.GetParameterSet(), resistor);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box->SetCoord(2 * n, reverse ? stop[n] : start[n]);
			box->SetCoord(2 * n + 1, reverse ? start[n] : stop[n]);
		}
		signal.SetupDiracPulse(30e9);
		if (compact)
		{
			const std::vector<double> region{fine[0][2], fine[1][2], fine[2][2], fine[0][6], fine[1][6], fine[2][6]};
			op.reset(Operator_CartesianSubgrid::New(region, refined ? Grid::Index{2, 2, 2} : Grid::Index{1, 1, 1}));
		}
		else
		{
			op.reset(Operator::New(1));
		}
		op->SetEngineType(Engine::BASIC);
		if (compact)
		{
			op->SetTimestepFactor(.5);
		}
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "rectangular resistor fixture geometry failed");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(bc);
	}
	Operator_CartesianSubgrid &Compact()
	{
		auto *value = dynamic_cast<Operator_CartesianSubgrid *>(op.get());
		Require(value, "expected compact operator");
		return *value;
	}
};
static double Seed(unsigned int n, const Grid::Point &p)
{
	const double pi = std::acos(-1.);
	return std::sin(pi * p[(n + 1) % 3] / .004) * std::sin(pi * p[(n + 2) % 3] / .004) *
		std::cos(2 * pi * p[n] / .004);
}
static void Control(unsigned int direction, bool graded, bool sheet, bool caps, bool reverse,
	const SubgridDeviceTarget *target, bool domains, bool refined)
{
	Fixture native(true, graded, direction, sheet, caps, reverse, refined), ordinary(false, graded, direction, sheet, caps, reverse, refined), reference(true, graded, direction, sheet, caps, reverse, refined);
	Require(native.op->CalcECOperator() == 0, "native rectangular resistor preparation failed");
	const double dt = native.op->GetTimestep();
	Require(reference.op->CalcECOperator() == 0 && reference.op->GetTimestep() == dt, "native CPU reference preparation differs");
	ordinary.op->SetTimestep(dt);
	Require(ordinary.op->CalcECOperator() == 0 && ordinary.op->GetTimestepValid(), "ordinary rectangular resistor preparation or shared timestep failed");
	auto &grid = native.Compact().Grid();
	Grid::Point a = native.start, b = native.stop;
	for (unsigned int n = 0; n < 3; ++n)
	{
		a[n] *= .001;
		b[n] *= .001;
	}
	const auto cells = grid.DistributedResistor(reverse ? b : a, reverse ? a : b, direction, 50);
	Require(cells.size() == (refined ? (sheet ? 6 : 18) : (sheet ? 42 : 294)), "native distributed edge count changed");
	double conductanceError = 0;
	for (const auto &cell : cells)
	{
		unsigned int pos[3];
		for (unsigned int n = 0; n < 3; ++n)
		{
			auto node = std::min_element(ordinary.fine[n].begin(), ordinary.fine[n].end(), [&](double x, double y)
			{
				return std::abs(x - cell.start[n] / .001) < std::abs(y - cell.start[n] / .001);
			});
			pos[n] = static_cast<unsigned int>(node - ordinary.fine[n].begin());
		}
		const double vv = ordinary.op->GetVV(direction, pos), vi = ordinary.op->GetVI(direction, pos);
		const double g = (1 - vv) / vi;
		conductanceError = std::max(conductanceError, std::abs(g - 1 / cell.resistance) * cell.resistance);
	}
	Require(conductanceError < 2e-5, "distributed conductance differs from actual ordinary Yee coefficients");
	if (target)
	{
		native.op->SetEngineType(Engine::GPU);
		if (domains)
		{
			std::vector<std::unique_ptr<GPU_Backend>> backends;
			for (unsigned int worker = 0; worker < 2; ++worker)
			{
				target->SelectDomain(worker);
				backends.emplace_back(GPU_Backend::New(target->backend));
				Require(bool(backends.back()), "rectangular-port domain backend unavailable");
				target->Verify(*backends.back());
			}
			native.Compact().SetDomainExecution((direction + 1) % 3, {1.6}, std::move(backends));
		}
		else
		{
			target->SelectDomain(0);
			std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target->backend));
			Require(bool(backend), "rectangular-port backend unavailable");
			target->Verify(*backend);
			native.op->SetGPUBackend(backend.release());
		}
	}
	grid.SetElectric(Seed);
	reference.Compact().Grid().SetElectric(Seed);
	reference.engine.reset(reference.op->CreateEngine());
	native.engine.reset(native.op->CreateEngine());
	ordinary.engine.reset(ordinary.op->CreateEngine());
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int x = 0; x < 9; ++x)
		{
			for (unsigned int y = 0; y < 9; ++y)
			{
				for (unsigned int z = 0; z < 9; ++z)
				{
					unsigned int pos[3] = {x, y, z};
					if (pos[n] == 8 || (ordinary.op->GetVV(n, pos) == 0 && ordinary.op->GetVI(n, pos) == 0))
					{
						continue;
					}
					double p[3];
					ordinary.op->GetYeeCoords(n, pos, p, false);
					ordinary.engine->SetVolt(n, pos, Seed(n, {p[0] * .001, p[1] * .001, p[2] * .001}) *
						ordinary.op->GetEdgeLength(n, pos));
				}
			}
		}
	}
	const double initial = grid.Energy(dt);
	Require(initial > 0, "distributed-port initial energy is empty");
	double peak = 0, error = 0, balance = 0;
	for (unsigned int batch = 0; batch < 8; ++batch)
	{
		native.engine->IterateTS(16);
		ordinary.engine->IterateTS(16);
		reference.engine->IterateTS(16);
		double loss = 0;
		for (std::size_t k = 0; k < native.engine->GetExtensionCount(); ++k)
		{
			if (auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(native.engine->GetExtension(k)))
			{
				loss += extension->ElectricLoss() + extension->MagneticLoss() + extension->SurfaceLoss();
			}
		}
		Require(loss > 0, "distributed resistor dissipates no power");
		balance = std::max(balance, std::abs(grid.Energy(dt) + loss - initial) / initial);
		for (unsigned int x = 1; x < 8; ++x)
		{
			for (unsigned int y = 1; y < 8; ++y)
			{
				for (unsigned int z = 1; z < 8; ++z)
				{
					unsigned int pos[3] = {x, y, z};
					for (unsigned int n = 0; n < 3; ++n)
					{
						double wanted = ordinary.engine->GetVolt(n, pos) / ordinary.op->GetEdgeLength(n, pos);
						double p[3];
						ordinary.op->GetYeeCoords(n, pos, p, false);
						const Grid::Point point{p[0] * .001, p[1] * .001, p[2] * .001};
						if (refined)
						{
							wanted = reference.Compact().Grid().Electric(n, point);
						}
						const double actual = grid.Electric(n, point);
						peak = std::max(peak, std::abs(wanted));
						error = std::max(error, std::abs(actual - wanted));
					}
				}
			}
		}
	}
	Require(peak > 0 && error / peak < 2e-4 && balance < (target ? 2e-4 : 2e-10),
		"distributed-port fields or passive energy differ from ordinary Yee");
	std::cout << "NATIVE_BOX_PORT axis=" << direction << " graded=" << graded << " sheet=" << sheet
		<< " refined=" << refined << " caps=" << caps << " reverse=" << reverse << " device=" << bool(target) << " domains=" << domains
		<< " field_reference=" << (refined ? "nativeCPU" : "ordinaryYee") << " edges=" << cells.size() << " conductance_error=" << conductanceError
		<< " field_error=" << error / peak << " balance=" << balance << std::endl;
	if (!caps)
	{
		a[direction] += .00001;
		Reject([&] { grid.DistributedResistor(a, b, direction, 50); });
	}
}
int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "expected optional backend flag");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool graded : {false, true})
			{
				for (bool sheet : {false, true})
				{
					const bool caps = sheet;
					for (bool refined : {false, true})
					{
						Control(axis, graded, sheet, caps, graded, nullptr, false, refined);
						Control(axis, graded, sheet, caps, graded, &target, target.physicalDevices > 1, refined);
					}
				}
			}
		}
		std::cout << "PASS: native rectangular R geometry, actual ordinary coefficients, passive fields and device domains" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

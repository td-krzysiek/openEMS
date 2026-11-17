/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "tools/constants.h"
#include "ContinuousStructure.h"
#include "CSPrimCylinder.h"
#include "CSPrimLinPoly.h"
#include "CSPropMaterial.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "subgrid_device_target.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
using Advance = Engine_Ext_CartesianSubgrid_Conductivity;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

// Analytic membership and an explicit midpoint sum deliberately do not call
// CSXCAD, its priority resolver or the implementation's quadrature helper.
static double Fraction(const Grid::Point &lo, const Grid::Point &hi, unsigned int axis)
{
	const unsigned int a = (axis + 1) % 3, b = (axis + 2) % 3;
	unsigned int inside = 0;
	for (unsigned int i = 0; i < 4; ++i)
	{
		for (unsigned int j = 0; j < 4; ++j)
		{
			for (unsigned int k = 0; k < 4; ++k)
			{
				const unsigned int indices[3] = {i, j, k};
				Grid::Point p;
				for (unsigned int n = 0; n < 3; ++n)
				{
					p[n] = 1000 * (lo[n] + (hi[n] - lo[n]) * (indices[n] + .5) / 4);
				}
				const double u = p[a], v = p[b], z = p[axis];
				const bool slab = z > 1 && z < 5 && u > 1 && u < 5 && v > 1 && v < 5 &&
					u + v > 2.6 && u + v < 9.4 && u - v < 3.4 && v - u < 3.4;
				const bool bore = z > 1 && z < 5 && std::hypot(u - 3, v - 3) < .6;
				inside += slab && !bore;
			}
		}
	}
	return double(inside) / 64;
}

struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	std::unique_ptr<Engine> engine;
	Advance *advance = nullptr;

	Fixture(bool native, unsigned int axis, bool reverse, unsigned int mode)
	{
		const double unit = reverse ? 1e-6 : .001;
		const Grid::Point origin{31, -23, 17};
		const unsigned int a = (axis + 1) % 3, b = (axis + 2) % 3;
		const auto draw = [&](unsigned int n, double x) { return (origin[n] + x) * .001 / unit; };
		csx.GetGrid()->SetDeltaUnit(unit);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= 6; ++i)
			{
				csx.GetGrid()->AddDiscLine(n, draw(n, i));
			}
		}
		const auto add = [&](bool air)
		{
			auto *material = new CSPropMaterial(csx.GetParameterSet());
			material->SetEpsilon(air ? 1 : 3.556);
			material->SetMue(air ? 1 : 1.7);
			material->SetKappa(!air && native && (mode & 1) ? .02 : 0);
			material->SetSigma(!air && native && (mode & 2) ? 100 : 0);
			csx.AddProperty(material);
			if (air)
			{
				auto *cylinder = new CSPrimCylinder(csx.GetParameterSet(), material);
				cylinder->SetPriority(9);
				cylinder->SetRadius(.6 * .001 / unit);
				for (unsigned int n = 0; n < 3; ++n)
				{
					cylinder->SetCoord(2 * n, draw(n, n == axis ? 1 : 3));
					cylinder->SetCoord(2 * n + 1, draw(n, n == axis ? 5 : 3));
				}
			}
			else
			{
				auto *polygon = new CSPrimLinPoly(csx.GetParameterSet(), material);
				polygon->SetPriority(1);
				polygon->SetNormDir(int(axis));
				polygon->SetElevation(draw(axis, reverse ? 5 : 1));
				polygon->SetLength((reverse ? -4 : 4) * .001 / unit);
				const double points[8][2] = {{1.6, 1}, {4.4, 1}, {5, 1.6}, {5, 4.4},
					{4.4, 5}, {1.6, 5}, {1, 4.4}, {1, 1.6}};
				for (unsigned int i = 0; i < 8; ++i)
				{
					const auto &point = points[reverse ? 7 - i : i];
					polygon->AddCoord(draw(a, point[0]));
					polygon->AddCoord(draw(b, point[1]));
				}
			}
		};
		add(reverse);
		add(!reverse);
		Require(csx.Update().empty(), "shaped-loss native CAD is invalid");
		std::vector<double> region;
		for (double endpoint : {2., 4.})
		{
			for (unsigned int n = 0; n < 3; ++n)
			{
				region.push_back(draw(n, endpoint));
			}
		}
		signal.SetupDiracPulse(30e9);
		op.reset(Operator_CartesianSubgrid::New(region));
		op->SetMaxLeaves(4000000);
		op->SetMaterialSamples(4);
		op->SetEngineType(Engine::BASIC);
		op->SetExcitationSignal(&signal);
		Require(op->SetGeometryCSX(&csx), "shaped-loss geometry rejected");
		op->SetTimestepFactor(.7);
		if (!native)
		{
			op->AddExtension(new Loss(op.get(), [=](const Grid::Point &lo, const Grid::Point &hi)
			{
				const double fraction = Fraction(lo, hi, axis);
				return Loss::Coefficients{mode & 1 ? .02 * fraction : 0, mode & 2 ? 100 * fraction : 0};
			}));
		}
		Require(op->CalcECOperator() == 0 && op->GetNumberOfExtentions() == 1,
			"shaped-loss native/manual preparation or provider count differs");
		unsigned int partial = 0, air = 0;
		for (std::size_t leaf = 0; leaf < op->Grid().CellCount(); ++leaf)
		{
			const auto lo = op->Grid().LeafStart(leaf), length = op->Grid().LeafLengths(leaf);
			Grid::Point hi;
			for (unsigned int n = 0; n < 3; ++n)
			{
				hi[n] = lo[n] + length[n];
			}
			const double fraction = Fraction(lo, hi, axis);
			const auto medium = op->Grid().LeafMedium(leaf);
			Require(std::abs(medium.epsilon / (EPS0 * (1 + 2.556 * fraction)) - 1) < 2e-13 &&
				std::abs(medium.permeability / (MUE0 * (1 + .7 * fraction)) - 1) < 2e-13,
				"independent analytic material fraction disagrees before loss execution");
			partial += fraction > 0 && fraction < 1;
			air += fraction == 0;
		}
		Require(partial > 0 && air > 0, "shaped-loss test missed partial leaves or priority air");
		op->Grid().SetElectric([](unsigned int n, const Grid::Point &p)
		{
			return std::sin(1731 * p[0] + 933 * p[1] - 671 * p[2] + n);
		});
		op->Grid().SetMagnetic([](unsigned int n, const Grid::Point &p)
		{
			return std::cos(633 * p[0] - 731 * p[1] + 411 * p[2] + n) / 376.730313668;
		});
		engine.reset(op->CreateEngine());
		advance = dynamic_cast<Advance *>(engine->GetExtension(0));
		Require(advance, "shaped-loss paired extension missing");
	}
};

static double FieldError(const Grid &a, const Grid &b)
{
	const auto x = a.DeviceStateData(), y = b.DeviceStateData();
	double error = 0;
	for (const auto &pair : {std::make_pair(&x.electric, &y.electric), std::make_pair(&x.magnetic, &y.magnetic)})
	{
		Require(pair.first->size() == pair.second->size(), "shaped-loss field counts differ");
		double scale = 0, delta = 0;
		for (std::size_t i = 0; i < pair.first->size(); ++i)
		{
			Require(std::isfinite(pair.first->at(i)) && std::isfinite(pair.second->at(i)), "nonfinite shaped-loss field");
			scale = std::max(scale, std::abs(double(pair.first->at(i))));
			delta = std::max(delta, std::abs(double(pair.first->at(i)) - pair.second->at(i)));
		}
		error = std::max(error, delta / std::max(scale, 1e-300));
	}
	return error;
}

static void Control(unsigned int axis, bool reverse, unsigned int mode,
	const SubgridDeviceTarget &target, bool domains)
{
	Fixture cpu(true, axis, reverse, mode), manual(false, axis, reverse, mode),
		device(true, axis, reverse, mode), fp32(true, axis, reverse, mode);
	const double dt = cpu.op->GetTimestep(), initial = cpu.op->Grid().Energy(dt);
	Require(dt == manual.op->GetTimestep() && dt == device.op->GetTimestep(), "shaped-loss physical timestep changed");
	std::vector<std::unique_ptr<GPU_Backend>> owners;
	std::vector<GPU_Backend *> backends;
	for (std::size_t worker = 0; worker < (domains ? 3 : 1); ++worker)
	{
		target.SelectDomain(worker);
		owners.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(owners.back()), "shaped-loss device unavailable");
		target.Verify(*owners.back());
		std::cout << "SHAPED_LOSS_BACKEND owner=" << worker << " backend=" << owners.back()->GetName() << '\n';
		backends.push_back(owners.back().get());
	}
	std::unique_ptr<GPU_CartesianSubgrid> execution;
	std::shared_ptr<CartesianSubgridDomainPlan> plan;
	if (domains)
	{
		plan = std::make_shared<CartesianSubgridDomainPlan>(*device.advance, axis,
			std::vector<double>{.002, .004});
		execution = CreateCartesianSubgridDomains(plan, backends);
	}
	else
	{
		execution = device.advance->CreateDevice(*backends.front());
	}
	Require(bool(execution), "shaped-loss execution missing");
	std::unique_ptr<GPU_Backend> referenceBackend(GPU_Backend::New("reference"));
	Require(bool(referenceBackend) && referenceBackend->GetName().find("reference (") == 0,
		"independent shaped-loss FP32 reference unavailable");
	auto reference = fp32.advance->CreateDevice(*referenceBackend);
	Require(bool(reference), "shaped-loss FP32 execution missing");
	double cpuBalance = 0;
	for (unsigned int step = 0; step < 512; ++step)
	{
		cpu.advance->UpdateElectric(dt, 0);
		cpu.advance->UpdateMagnetic(dt);
		manual.advance->UpdateElectric(dt, 0);
		manual.advance->UpdateMagnetic(dt);
		execution->Step(dt, 0);
		reference->Step(dt, 0);
		if (step % 64 == 63)
		{
			cpuBalance = std::max(cpuBalance, std::abs(cpu.op->Grid().Energy(dt) +
				cpu.advance->ElectricLoss() + cpu.advance->MagneticLoss() - initial) / initial);
		}
	}
	execution->ReadFields(device.op->Grid());
	reference->ReadFields(fp32.op->Grid());
	fp32.advance->ReadDeviceDiagnostics(*reference);
	device.advance->ReadDeviceDiagnostics(*execution);
	const double manualError = FieldError(cpu.op->Grid(), manual.op->Grid());
	const double deviceError = FieldError(cpu.op->Grid(), device.op->Grid());
	const double fp32Error = FieldError(fp32.op->Grid(), device.op->Grid());
	const double fp32Ledger = std::max(
		std::abs(fp32.advance->ElectricLoss() - device.advance->ElectricLoss()),
		std::abs(fp32.advance->MagneticLoss() - device.advance->MagneticLoss())) / initial;
	Require(fp32Error < 2e-6 && fp32Ledger < 2e-6,
		"shaped-loss physical/reference FP32 fields or ledger exceeded retained limits");
	double manualLedger = 0, deviceLedger = 0;
	for (const auto &pair : {std::make_pair(cpu.advance->ElectricLoss(), manual.advance->ElectricLoss()),
		std::make_pair(cpu.advance->MagneticLoss(), manual.advance->MagneticLoss())})
	{
		manualLedger = std::max(manualLedger, std::abs(pair.first - pair.second) / initial);
	}
	for (const auto &pair : {std::make_pair(cpu.advance->ElectricLoss(), device.advance->ElectricLoss()),
		std::make_pair(cpu.advance->MagneticLoss(), device.advance->MagneticLoss())})
	{
		deviceLedger = std::max(deviceLedger, std::abs(pair.first - pair.second) / initial);
	}
	const double balance = std::abs(device.op->Grid().Energy(dt) + device.advance->ElectricLoss() +
		device.advance->MagneticLoss() - initial) / initial;
	Require(manualError < 2e-8 && manualLedger < 2e-10 && cpuBalance < 2e-10,
		"native shaped-loss mapping differs from independent analytic loss or CPU passivity");
	Require(deviceError < 2e-4 && deviceLedger < 2e-4 && balance < 2e-4,
		"shaped-loss device field/ledger/passivity exceeded retained limits");
	Require((mode & 1 ? device.advance->ElectricLoss() > 0 : device.advance->ElectricLoss() == 0) &&
		(mode & 2 ? device.advance->MagneticLoss() > 0 : device.advance->MagneticLoss() == 0),
		"shaped-loss conductivity mode was dropped or spuriously added");
	std::cout << "SHAPED_LOSS axis=" << axis << " reversed=" << reverse << " mode=" << mode
		<< " workers=" << backends.size() << " steps=512 manual_error=" << manualError
		<< " manual_ledger=" << manualLedger << " cpu_balance=" << cpuBalance
		<< " device_error=" << deviceError << " device_ledger=" << deviceLedger
		<< " balance=" << balance << " fp32_error=" << fp32Error << " fp32_ledger=" << fp32Ledger
		<< " electric_loss=" << device.advance->ElectricLoss() / initial
		<< " magnetic_loss=" << device.advance->MagneticLoss() / initial << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		std::cout << std::setprecision(17);
		Require(argc <= 2, "unexpected shaped-loss arguments");
		const bool referenceDomains = argc == 2 && std::string(argv[1]) == "--reference-domains";
		const auto target = argc == 1 || referenceDomains ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
		const bool domains = referenceDomains || target.physicalDevices > 1;
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool reverse : {false, true})
			{
				for (unsigned int mode = 1; mode <= 3; ++mode)
				{
					Control(axis, reverse, mode, target, domains);
				}
			}
		}
		std::cout << "PASS: native shaped electric/magnetic conductivity, analytic priority fractions and passive device domains\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

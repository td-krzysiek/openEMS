/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "CSPropMaterial.h"
#include "CSPrimBox.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "subgrid_device_target.h"
#include "tools/constants.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using Loss = Operator_Ext_CartesianSubgrid_Conductivity;
static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}
struct Configuration
{
	unsigned int normal = 0;
	int direction = 1;
	double electric = .05, magnetic = 0;
	bool native = false, lossFirst = false, addLoss = true, overlap = false, guardOnly = false;
	bool soft = false, ports = false, corner = false, surface = false;
};
struct Fixture
{
	ContinuousStructure csx;
	Excitation signal;
	std::unique_ptr<Operator_CartesianSubgrid> op;
	Fixture(const Configuration &test)
	{
		csx.GetGrid()->SetDeltaUnit(.001);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i <= (test.corner ? 24u : n == test.normal ? 80u : 4u); ++i)
			{
				csx.GetGrid()->AddDiscLine(n, i);
			}
		}
		signal.SetupDiracPulse(30e9);
		std::vector<double> box{1, 1, 1, 3, 3, 3};
		box[test.normal] = 24;
		box[test.normal + 3] = 40;
		CartesianSubgrid::Index ratio{1, 1, 1};
		ratio[test.normal] = 2;
		if (test.corner)
		{
			box = {10, 10, 10, 14, 14, 14};
			ratio = {2, 2, 2};
		}
		op.reset(Operator_CartesianSubgrid::New(box, ratio));
		op->SetExcitationSignal(&signal);
		op->SetEngineType(Engine::BASIC);
		Require(op->SetGeometryCSX(&csx), "joint PML geometry rejected");
		int bc[6] = {0, 0, 0, 0, 0, 0};
		const auto hn = (test.normal + 2) % 3;
		bc[2 * hn] = bc[2 * hn + 1] = 1;
		bc[2 * test.normal] = bc[2 * test.normal + 1] = 3;
		if (test.corner)
		{
			std::fill(bc, bc + 6, 3);
		}
		op->SetBoundaryCondition(bc);
		op->SetTimestepFactor(.8);
		auto sample = [test](const CartesianSubgrid::Point &lo, const CartesianSubgrid::Point &hi)
		{
			const double position = .5 * (lo[test.normal] + hi[test.normal]);
			bool inside = test.overlap || (test.guardOnly ? position >= .008 && position < .009 :
				position >= .032 && position < .040);
			if (test.corner)
			{
				inside = true;
				for (unsigned int n = 0; n < 3; ++n)
				{
					const double coordinate = .5 * (lo[n] + hi[n]);
					inside = inside && coordinate >= .010 && coordinate < .014;
				}
			}
			return inside ? Loss::Coefficients{test.electric, test.magnetic} : Loss::Coefficients{};
		};
		if (test.native)
		{
			auto *material = new CSPropMaterial(csx.GetParameterSet());
			material->SetKappa(test.electric);
			material->SetSigma(test.magnetic);
			csx.AddProperty(material);
			auto *materialBox = new CSPrimBox(csx.GetParameterSet(), material);
			for (unsigned int n = 0; n < 3; ++n)
			{
				materialBox->SetCoord(2 * n, test.corner ? 10. : n == test.normal && !test.overlap ? 32. : 0.);
				materialBox->SetCoord(2 * n + 1, test.corner ? 14. : n == test.normal ? (test.overlap ? 80. : 40.) : 4.);
			}
		}
		auto *loss = test.addLoss && !test.native ? new Loss(op.get(), sample) : nullptr;
		if (loss && test.surface)
		{
			const auto en = (test.normal + 1) % 3;
			Loss::SurfaceAdmittance surface;
			surface.start = {.002, .002, .002};
			surface.start[test.normal] = test.corner ? .012 : .036;
			surface.start[en] = .001;
			surface.stop = surface.start;
			surface.stop[en] = .002;
			if (test.corner)
			{
				surface.start = {.012, .012, .012};
				surface.start[(en + 1) % 3] += .0005;
				surface.stop = surface.start;
				surface.stop[en] += .0005;
			}
			surface.conductance = .01;
			surface.poles = {{25, 1e-10}, {50, 2e-10}};
			loss->SetSurfaceAdmittances({surface});
		}
		if (loss && test.lossFirst)
		{
			op->AddExtension(loss);
		}
		std::array<unsigned int, 6> cells{};
		cells[2 * test.normal] = cells[2 * test.normal + 1] = 8;
		if (test.corner)
		{
			cells.fill(4);
		}
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), cells));
		if (loss && !test.lossFirst)
		{
			op->AddExtension(loss);
		}
		Require(op->CalcECOperator() == 0, "joint PML preparation failed");
		for (std::size_t i = 0; i < op->GetNumberOfExtentions(); ++i)
		{
			if (const auto *pml = dynamic_cast<Operator_Ext_CartesianSubgrid_PML *>(op->GetExtension(i)))
			{
				const auto &footprint = pml->PreparedFootprint();
				Require(!footprint.electric.empty() && !footprint.magnetic.empty(), "empty protected PML support");
			}
		}
	}
};
struct Result
{
	std::vector<double> electric, magnetic;
	double balance = 0, loss = 0, residual = 0;
	double sourceWork = 0, surfaceLoss = 0;
};
static Result Run(const Configuration &test, const SubgridDeviceTarget *device)
{
	Fixture fixture(test);
	auto &op = *fixture.op;
	const double dt = op.GetTimestep(), centre = test.direction > 0 ? .020 : .060;
	const unsigned int en = (test.normal + 1) % 3, hn = (test.normal + 2) % 3;
	const double impedance = std::sqrt(MUE0 / EPS0);
	auto pulse = [&](double position) { return std::exp(-std::pow((position - centre) / .004, 2)); };
	op.Grid().SetElectric([&](unsigned int n, const CartesianSubgrid::Point &p)
	{
		if (test.corner)
		{
			const double x = (p[0] - .012) / .002, y = (p[1] - .012) / .002, z = (p[2] - .012) / .002;
			return (n == 1 ? -z : n == 2 ? y : 0) * std::exp(-x * x - y * y - z * z);
		}
		return n == en ? pulse(p[test.normal]) : 0;
	});
	op.Grid().SetMagnetic([&](unsigned int n, const CartesianSubgrid::Point &p)
	{
		if (test.corner)
		{
			return 0.;
		}
		return n == hn ? test.direction * pulse(p[test.normal] + test.direction * C0 * dt / 2) / impedance : 0;
	});
	std::size_t selected = op.Grid().EdgeCount();
	if (test.soft)
	{
		CartesianSubgrid::Point position{.002, .0015, .002};
		position[test.normal] = .036;
		double distance = 1;
		for (std::size_t edge = 0; edge < op.Grid().EdgeCount(); ++edge)
		{
			if (op.Grid().EdgeDirection(edge) != en)
			{
				continue;
			}
			const auto point = op.Grid().EdgePosition(edge);
			double square = 0;
			for (unsigned int n = 0; n < 3; ++n)
			{
				square += std::pow(point[n] - position[n], 2);
			}
			if (square < distance)
			{
				distance = square;
				selected = edge;
			}
		}
		Require(selected < op.Grid().EdgeCount(), "joint source edge missing");
	}
	if (test.ports)
	{
		CartesianSubgrid::Point a{.002, .002, .002}, b = a, c = a, d = a;
		a[test.normal] = b[test.normal] = .036;
		c[test.normal] = d[test.normal] = .037;
		a[en] = c[en] = .001;
		b[en] = d[en] = .002;
		op.Grid().SetResistors({{a, b, 50, 1e-5}, {c, d, 75, 0}}, dt);
		for (std::size_t i = 0; i < op.GetNumberOfExtentions(); ++i)
		{
			if (auto *loss = dynamic_cast<Loss *>(op.GetExtension(i)))
			{
				Require(loss->BuildExtension(), "joint port responses did not rebuild");
				loss->ValidatePreparedExtensions();
			}
		}
	}
	if (test.soft)
	{
		// Ports must be configured before soft-source metadata. SetSources
		// accepts physical E, then applies sqrt(mass) internally.
		op.Grid().SetSources({{selected, .5}});
	}
	if (test.soft || test.ports)
	{
		fixture.signal.SetupDiracPulse(30e9);
		Require(fixture.signal.buildExcitationSignal(1024) && fixture.signal.GetLength() > 0 &&
			std::any_of(fixture.signal.GetVoltageSignal(),
				fixture.signal.GetVoltageSignal() + fixture.signal.GetLength(), [](double value) { return value != 0; }),
			"joint test drive is empty");
	}
	if (device)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(device->backend));
		Require(bool(backend), "joint PML backend unavailable");
		device->Verify(*backend);
		op.SetEngineType(Engine::GPU);
		op.SetGPUBackend(backend.release());
	}
	std::unique_ptr<Engine> engine(op.CreateEngine());
	const double initial = op.Grid().Energy(dt);
	const auto steps = test.corner ? 1000u : static_cast<unsigned int>(std::ceil(5 * .080 / (C0 * dt)));
	Result result;
	for (unsigned int step = 0; step < steps; ++step)
	{
		engine->IterateTS(1);
		double boundary = 0, material = 0, history = 0;
		for (std::size_t i = 0; i < engine->GetExtensionCount(); ++i)
		{
			if (const auto *pml = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(engine->GetExtension(i)))
			{
				boundary += pml->BoundaryWorkOut();
			}
			if (const auto *loss = dynamic_cast<Engine_Ext_CartesianSubgrid_Conductivity *>(engine->GetExtension(i)))
			{
				material += loss->ElectricLoss() + loss->MagneticLoss() + loss->SurfaceLoss();
				history += loss->AdditionalStoredEnergy();
				result.surfaceLoss = loss->SurfaceLoss() / initial;
			}
		}
		const double energy = op.Grid().Energy(dt) + history;
		const double available = initial + std::max(0., op.Grid().SourceWork());
		Require(std::isfinite(energy) && material >= 0 && energy / available < 1.1,
			"joint PML evolution exceeds initial and supplied energy");
		result.balance = std::max(result.balance,
			std::abs(energy + material + boundary + op.Grid().ResistorLoss() - op.Grid().SourceWork() - initial) / initial);
		result.sourceWork = op.Grid().SourceWork() / initial;
		result.loss = material / initial;
		result.residual = energy / initial;
		if (step % 16 == 0)
		{
			for (unsigned int sample = 0; sample < 16; ++sample)
			{
				CartesianSubgrid::Point point{.002, .002, .002};
				point[test.normal] = .002 + sample * .005;
				if (test.corner)
				{
					for (unsigned int n = 0; n < 3; ++n)
					{
						point[n] = sample & (1u << n) ? (sample >= 8 ? .022 : .018) :
							(sample >= 8 ? .002 : .006);
					}
				}
				for (unsigned int n = 0; n < (test.corner ? 3u : 1u); ++n)
				{
					result.electric.push_back(op.Grid().Electric(test.corner ? n : en, point));
					result.magnetic.push_back(op.Grid().Magnetic(test.corner ? n : hn, point));
				}
			}
		}
	}
	Require(result.balance < (device ? 2e-4 : 2e-9), "joint PML/material power identity failed");
	Require(result.residual < (test.corner ? .05 : .01), "joint PML did not absorb the outgoing packet");
	if (test.addLoss && (test.electric != 0 || test.magnetic != 0))
	{
		Require(result.loss > 1e-4, "joint packet did not interact with material loss");
	}
	if (test.soft || test.ports)
	{
		Require(std::abs(result.sourceWork) > 1e-8, "joint drive did no measured work");
	}
	if (test.surface)
	{
		Require(result.surfaceLoss > 1e-5, "joint packet did not interact with surface histories");
	}
	std::cout << "JOINT_PML normal=" << test.normal << " direction=" << test.direction << " native=" << test.native
		<< " first=" << test.lossFirst << " device=" << bool(device) << " balance=" << result.balance
		<< " material_loss=" << result.loss << " residual=" << result.residual << " soft=" << test.soft
		<< " ports=" << test.ports << " corner=" << test.corner << " source_work=" << result.sourceWork
		<< " surface=" << test.surface << " surface_loss=" << result.surfaceLoss << std::endl;
	return result;
}
static double Compare(const Result &a, const Result &b, bool exact = false)
{
	double error = 0;
	for (const auto &pair : {std::make_pair(&a.electric, &b.electric), std::make_pair(&a.magnetic, &b.magnetic)})
	{
		Require(pair.first->size() == pair.second->size(), "joint packet dimensions differ");
		double scale = 0, delta = 0;
		for (std::size_t i = 0; i < pair.first->size(); ++i)
		{
			scale = std::max(scale, std::abs(pair.first->at(i)));
			delta = std::max(delta, std::abs(pair.first->at(i) - pair.second->at(i)));
		}
		error = std::max(error, delta / std::max(scale, 1e-300));
	}
	Require(error <= (exact ? 0 : 2e-4), "joint PML packet fields differ");
	return error;
}
int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 3, "expected optional --surface-only and --metal/--cuda");
		SubgridDeviceTarget target;
		bool surfaceOnly = false;
		for (int argument = 1; argument < argc; ++argument)
		{
			if (std::string(argv[argument]) == "--surface-only")
			{
				surfaceOnly = true;
			}
			else
			{
				target = SubgridDeviceTarget::FromFlag(argv[argument]);
			}
		}
		if (surfaceOnly)
		{
			for (bool corner : {false, true})
			{
				Configuration test;
				test.surface = true;
				test.corner = corner;
				test.soft = !corner;
				test.ports = !corner;
				test.magnetic = test.electric * MUE0 / EPS0;
				const auto cpu = Run(test, nullptr);
				const double error = Compare(cpu, Run(test, &target));
				std::cout << "JOINT_SURFACE corner=" << corner << " field_error=" << error << std::endl;
			}
			std::cout << "PASS: disjoint surface-state/PML fields, stored energy and work" << std::endl;
			return 0;
		}
		for (unsigned int normal : {0u, 1u, 2u})
		{
			Configuration test;
			test.normal = normal;
			test.direction = normal == 1 ? -1 : 1;
			test.magnetic = test.electric * MUE0 / EPS0;
			const auto cpu = Run(test, nullptr);
			test.lossFirst = true;
			Compare(cpu, Run(test, nullptr), true);
			test.native = true;
			Compare(cpu, Run(test, nullptr), true);
			const double error = Compare(cpu, Run(test, &target));
			std::cout << "JOINT_COMPARE normal=" << normal << " field_error=" << error << std::endl;
		}
		Configuration zero;
		zero.electric = zero.magnetic = 0;
		const auto withZero = Run(zero, &target);
		zero.addLoss = false;
		Compare(withZero, Run(zero, &target), true);
		for (bool magnetic : {false, true})
		{
			Configuration bad;
			bad.electric = magnetic ? 0 : .02;
			bad.magnetic = magnetic ? 100 : 0;
			for (bool guardOnly : {false, true})
			{
				bad.overlap = !guardOnly;
				bad.guardOnly = guardOnly;
				bool rejected = false;
				try
				{
					Fixture fixture(bad);
				}
				catch (const std::invalid_argument &error)
				{
					rejected = std::string(error.what()).find("support overlap") != std::string::npos;
				}
				Require(rejected, "actual PML target/collar material overlap was not rejected");
			}
		}
		for (bool ports : {false, true})
		{
			Configuration driven;
			driven.magnetic = driven.electric * MUE0 / EPS0;
			driven.soft = true;
			driven.ports = ports;
			const auto cpu = Run(driven, nullptr);
			const double error = Compare(cpu, Run(driven, &target));
			std::cout << "JOINT_DRIVEN ports=" << ports << " field_error=" << error << std::endl;
		}
		for (bool magnetic : {false, true})
		{
			Configuration pure;
			pure.electric = magnetic ? 0 : .05;
			pure.magnetic = magnetic ? .05 * MUE0 / EPS0 : 0;
			const auto cpu = Run(pure, nullptr);
			const double error = Compare(cpu, Run(pure, &target));
			std::cout << "JOINT_PURE magnetic=" << magnetic << " field_error=" << error << std::endl;
		}
		Configuration corner;
		corner.corner = true;
		corner.native = true;
		corner.magnetic = corner.electric * MUE0 / EPS0;
		const auto cpuCorner = Run(corner, nullptr);
		const double cornerError = Compare(cpuCorner, Run(corner, &target));
		std::cout << "JOINT_CORNER field_error=" << cornerError << std::endl;
		Fixture ready(Configuration{});
		Operator_Ext_CartesianSubgrid_PML unprepared(ready.op.get(), {8, 8, 0, 0, 0, 0});
		bool rejectedUnprepared = false;
		try
		{
			Engine_Ext_CartesianSubgrid_PML invalid(&unprepared);
		}
		catch (const std::invalid_argument &error)
		{
			rejectedUnprepared = std::string(error.what()).find("before successful preparation") != std::string::npos;
		}
		Require(rejectedUnprepared, "unprepared PML engine was admitted");
		std::cout << "PASS: disjoint native/material PML, priority order, zero-loss exactness and overlap guards" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

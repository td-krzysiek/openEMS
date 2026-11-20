/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

static void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::runtime_error(message);
	}
}

struct Case
{
	unsigned int normal = 0;
	int direction = 1;
	bool refined = false, graded = false, absorb = true, corner = false, multiple = false;
	unsigned int minimumSteps = 0;
};

struct Result
{
	std::vector<double> electric, magnetic;
	double ledger = 0, returned = 0, residual = 0, peak = 0, timestep = 0;
};

static Result Run(const Case &test, const SubgridDeviceTarget *device)
{
	const unsigned int en = (test.normal + 1) % 3, hn = (test.normal + 2) % 3;
	const double speed = 1 / std::sqrt(epsilon * mu), impedance = std::sqrt(mu / epsilon);
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
	{
		const unsigned int count = test.corner ? 24 : (n == test.normal ? 80 : 4);
		for (unsigned int i = 0; i <= count; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, i - (test.graded && n == test.normal && i % 2 ? .15 : 0));
		}
	}
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op;
	if (test.corner)
	{
		op.reset(test.multiple ? Operator_CartesianSubgrid::NewRegions(
			{{8, 10, 10, 10, 14, 14}, {14, 10, 10, 16, 14, 14}}, {{2, 1, 2}, {1, 2, 2}}) :
			Operator_CartesianSubgrid::New({10, 10, 10, 14, 14, 14}));
	}
	else
	{
		std::vector<double> box{1, 1, 1, 3, 3, 3};
		box[test.normal] = 24;
		box[test.normal + 3] = 40;
		Grid::Index ratio{1, 1, 1};
		ratio[test.normal] = test.refined ? 2 : 1;
		op.reset(Operator_CartesianSubgrid::New(box, ratio));
	}
	if (device)
	{
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(device->backend));
		Require(bool(backend), "requested PML test backend unavailable");
		device->Verify(*backend);
		std::cout << "PML backend: " << backend->GetName() << std::endl;
		op->SetGPUBackend(backend.release());
		op->SetEngineType(Engine::GPU);
	}
	else
	{
		op->SetEngineType(Engine::BASIC);
	}
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "PML test geometry failed");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	std::array<unsigned int, 6> cells{};
	if (test.corner)
	{
		std::fill(bc, bc + 6, 3);
		cells.fill(4);
	}
	else
	{
		bc[2 * hn] = bc[2 * hn + 1] = 1;
		if (test.absorb)
		{
			bc[2 * test.normal] = bc[2 * test.normal + 1] = 3;
			cells[2 * test.normal] = cells[2 * test.normal + 1] = 8;
		}
	}
	op->SetBoundaryCondition(bc);
	if (test.absorb)
	{
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), cells));
	}
	Require(op->CalcECOperator() == 0, "PML test coefficients failed");
	const double dt = (test.corner ? 1 : .8) * op->GetTimestep();
	op->SetTimestep(dt);
	Require(op->CalcECOperator() == 0, "PML test timestep rebuild failed");
	if (test.corner)
	{
		op->Grid().SetElectric([](unsigned int n, const Grid::Point &p)
		{
			const double x = (p[0] - .012) / .002, y = (p[1] - .012) / .002, z = (p[2] - .012) / .002;
			return (n == 1 ? -z : n == 2 ? y : 0) * std::exp(-x * x - y * y - z * z);
		});
	}
	else
	{
		const double centre = test.direction > 0 ? .020 : .060;
		const auto pulse = [&](double position) { return std::exp(-std::pow((position - centre) / .004, 2)); };
		op->Grid().SetElectric([&](unsigned int n, const Grid::Point &p)
		{
			return n == en ? pulse(p[test.normal]) : 0;
		});
		op->Grid().SetMagnetic([&](unsigned int n, const Grid::Point &p)
		{
			return n == hn ? test.direction * pulse(p[test.normal] + test.direction * speed * dt / 2) / impedance : 0;
		});
	}
	std::unique_ptr<Engine> engine(op->CreateEngine());
	Result result;
	result.timestep = dt;
	const double initial = op->Grid().Energy(dt);
	result.peak = 1;
	Grid::Point probe{.002, .002, .002};
	probe[test.normal] = test.direction > 0 ? .016 : .064;
	const unsigned int steps = std::max(test.minimumSteps,
		test.corner ? 1000u : unsigned(5 * .080 / (speed * dt)));
	for (unsigned int step = 0; step < steps; ++step)
	{
		engine->IterateTS(1);
		double work = 0;
		for (std::size_t i = 0; i < engine->GetExtensionCount(); ++i)
		{
			const auto *paired = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(engine->GetExtension(i));
			Require(paired != nullptr, "unexpected PML test extension");
			work += paired->BoundaryWorkOut();
		}
		const double energy = op->Grid().Energy(dt);
		result.ledger = std::max(result.ledger, std::abs(energy + work - initial) / initial);
		result.peak = std::max(result.peak, energy / initial);
		if (!test.corner && (step + 1) * dt > 1.3 * .080 / speed)
		{
			result.returned = std::max(result.returned, std::abs(op->Grid().Electric(en, probe)));
		}
		if ((step + 1) % 8 == 0 || step + 1 == steps)
		{
			for (unsigned int position = 0; position < 16; ++position)
			{
				Grid::Point point{.002, .002, .002};
				if (test.corner)
				{
					for (unsigned int n = 0; n < 3; ++n)
					{
						point[n] = position & (1u << n) ? (position >= 8 ? .022 : .018) :
							(position >= 8 ? .002 : .006);
					}
				}
				else
				{
					point[test.normal] = .002 + position * .005;
				}
				for (unsigned int n = 0; n < 3; ++n)
				{
					result.electric.push_back(op->Grid().Electric(n, point));
					result.magnetic.push_back(op->Grid().Magnetic(n, point));
				}
			}
		}
	}
	result.residual = op->Grid().Energy(dt) / initial;
	Require(result.ledger < (device ? 2e-4 : 1e-9), "PML signed work balance failed");
	Require(std::isfinite(result.residual) && result.peak < 1.1, "PML test is unbounded");
	if (test.corner)
	{
		Require(result.residual < .05, "PML corner failed to absorb");
	}
	else if (test.absorb)
	{
		Require(result.returned < .03 && result.residual < .005, "PML packet reflection/decay failed");
	}
	else
	{
		Require(result.returned > .5 && result.residual > .99, "PEC negative control failed");
	}
	return result;
}

static double Compare(const std::vector<double> &reference, const std::vector<double> &candidate)
{
	Require(reference.size() == candidate.size(), "PML field record dimensions differ");
	double scale = 0, delta = 0;
	for (std::size_t i = 0; i < reference.size(); ++i)
	{
		Require(std::isfinite(candidate[i]), "nonfinite PML field record");
		scale = std::max(scale, std::abs(reference[i]));
		delta = std::max(delta, std::abs(reference[i] - candidate[i]));
	}
	const double error = delta / std::max(scale, 1e-300);
	Require(error < 2e-4, "PML CPU/device field difference exceeds 2e-4");
	return error;
}

template <class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const std::invalid_argument &)
	{
		return;
	}
	throw std::runtime_error("invalid device PML/half-step state was accepted");
}

static void Guards(const SubgridDeviceTarget &target)
{
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "requested PML guard backend unavailable");
	target.Verify(*backend);
	Grid grid({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu);
	auto core = backend->CreateCartesianSubgrid(grid);
	const double dt = .8 * grid.StableTimestep();
	core->RequireFixedTimestep(dt);
	CartesianSubgridPMLData data;
	data.timestep = dt;
	data.electric.start = {0, 1};
	data.electric.source = {0};
	data.electric.target = {0};
	data.electric.group = {0, 1};
	data.electric.weight = {1};
	data.electric.low = {0};
	data.electric.decay = {.9f};
	data.electric.change = {-.1f};
	data.magnetic = data.electric;
	for (unsigned int kind = 0; kind < 5; ++kind)
	{
		auto invalid = data;
		if (kind == 0)
		{
			invalid.electric.group = {0, 0, 1};
		}
		else if (kind == 1)
		{
			invalid.electric.source[0] = static_cast<std::uint32_t>(core->Fields().magneticCount);
		}
		else if (kind == 2)
		{
			invalid.electric.decay[0] = 1.01f;
		}
		else if (kind == 3)
		{
			invalid.electric.low[0] = 1e-40f;
		}
		else
		{
			invalid.electric.start[1] = 2;
		}
		Reject([&] { GPU_CreateCartesianSubgridPML(*backend, core->Fields(), invalid); });
	}
	auto wrongStep = data;
	wrongStep.timestep *= .5;
	Reject([&] { GPU_CreateCartesianSubgridPML(*backend, core->Fields(), wrongStep); });
	Reject([&] { core->RequireFixedTimestep(.5 * dt); });
	auto pml = GPU_CreateCartesianSubgridPML(*backend, core->Fields(), data);
	Engine_Box whole{{0, 0, 0}, {8, 8, 8}};
	Reject([&] { pml->Hook(false, GPU_POST, whole, 0); });
	Reject([&] { core->UpdateMagnetic(dt); });
	pml->Hook(false, GPU_PRE, whole, 0);
	Reject([&] { pml->ReadBoundaryWorkOut(); });
	Reject([&] { pml->Hook(false, GPU_PRE, whole, 0); });
	core->UpdateElectric(dt, 0);
	Reject([&] { core->UpdateElectric(dt, 0); });
	Reject([&] { core->ReadFields(grid); });
	Reject([&] { core->UpdateMagnetic(.5 * dt); });
	pml->Hook(false, GPU_POST, whole, 0);
	pml->Hook(false, GPU_APPLY, whole, 0);
	pml->Hook(true, GPU_PRE, whole, 0);
	core->UpdateMagnetic(dt);
	pml->Hook(true, GPU_POST, whole, 0);
	pml->Hook(true, GPU_APPLY, whole, 0);
	core->ReadFields(grid);
	Require(pml->ReadBoundaryWorkOut() == 0, "unexcited PML guard injected work");
	Reject([&] { pml->Hook(false, GPU_PRE, whole, 2); });
	Reject([&] { core->RequireFixedTimestep(dt); });
	Reject([&] { core->Step(.5 * dt, 0); });
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		SubgridDeviceTarget target;
		bool longPackets = false, longCorners = false;
		for (int i = 1; i < argc; ++i)
		{
			const std::string flag = argv[i];
			if (flag == "--long")
			{
				longPackets = true;
			}
			else if (flag == "--long-corners")
			{
				longCorners = true;
			}
			else
			{
				target = SubgridDeviceTarget::FromFlag(flag);
			}
		}
		Guards(target);
		std::vector<Case> cases;
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (int sign : {-1, 1})
			{
				Case packet;
				packet.normal = n;
				packet.direction = sign;
				cases.push_back(packet);
			}
		}
		Case local;
		local.refined = true;
		local.minimumSteps = longPackets ? 100000 : 0;
		cases.push_back(local);
		local.graded = true;
		local.direction = -1;
		local.minimumSteps = 0;
		cases.push_back(local);
		Case negative;
		negative.absorb = false;
		cases.push_back(negative);
		Case corner;
		corner.corner = true;
		corner.minimumSteps = longCorners ? 10000 : 0;
		cases.push_back(corner);
		corner.multiple = true;
		corner.minimumSteps = 0;
		cases.push_back(corner);
		for (const auto &test : cases)
		{
			const auto cpu = Run(test, nullptr), device = Run(test, &target);
			Require(cpu.timestep == device.timestep, "PML CPU/device timesteps differ");
			const double e = Compare(cpu.electric, device.electric), h = Compare(cpu.magnetic, device.magnetic);
			std::cout << "PML device comparison normal=" << test.normal << " direction=" << test.direction
				<< " local=" << test.refined << " graded=" << test.graded << " corner=" << test.corner
				<< " multiple=" << test.multiple << " E_error=" << e << " H_error=" << h
				<< " ledger=" << device.ledger << " residual=" << device.residual << std::endl;
		}
		std::cout << "PASS: paired device PML fields, signed work, all normals/corners and PEC control" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

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
	throw std::runtime_error("invalid port-domain operation was accepted");
}

static Grid Fixture(unsigned int axis, bool drive, bool soft)
{
	Grid grid({12, 12, 12}, {.001, .001, .001}, {{3, 3, 3}, {9, 9, 9}}, epsilon, mu);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return .1 * std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .1 * std::cos(1381 * p[0] - 2171 * p[1] + 311 * p[2] + n) / 376.730313668;
	});
	grid.SetResistors({{{.001, .004, .004}, {.001, .004, .005}, 50, drive ? .001 : 0},
		{{.006, .006, .006}, {.006, .006, .0065}, 75, 0}}, .8 * grid.StableTimestep());
	if (soft)
	{
		const auto data = grid.DeviceData();
		const std::unordered_set<std::uint32_t> coupled(data.blockEdge.begin(), data.blockEdge.end());
		const CartesianSubgridDomainPlan plan(grid, axis, {.0047, .0083});
		bool found = false;
		for (const auto &domain : plan.Domains())
		{
			if (!domain.fields.rootConductance.empty())
			{
				continue;
			}
			for (std::size_t i = 0; i < domain.ownedElectric; ++i)
			{
				const auto edge = domain.globalElectric[i];
				if (coupled.count(edge))
				{
					grid.SetSources({{edge, 10}, {edge, -3}});
					found = true;
					break;
				}
			}
			if (found)
			{
				break;
			}
		}
		Require(found && !grid.DeviceData().sourceOffEdge.empty(),
			"no mass-coupled soft source on a worker without ports");
	}
	return grid;
}

static double FieldError(const Grid &a, const Grid &b)
{
	const auto first = a.DeviceData(), second = b.DeviceData();
	double error = 0;
	for (const auto &fields : {std::make_pair(&first.electric, &second.electric),
		std::make_pair(&first.magnetic, &second.magnetic)})
	{
		double scale = 0, delta = 0;
		Require(fields.first->size() == fields.second->size(), "field dimensions differ");
		for (std::size_t i = 0; i < fields.first->size(); ++i)
		{
			scale = std::max(scale, std::abs(double(fields.first->at(i))));
			delta = std::max(delta, std::abs(double(fields.first->at(i)) - fields.second->at(i)));
		}
		error = std::max(error, delta / std::max(scale, 1e-300));
	}
	return error;
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool drive,
	bool soft, bool mixed, bool dc, unsigned int steps)
{
	auto cpu = Fixture(axis, drive, soft);
	const double dt = .8 * cpu.StableTimestep();
	// Qualify attribution of a nonzero existing ledger, not just fresh zeros.
	for (unsigned int i = 0; i < 13; ++i)
	{
		cpu.Step(dt, .1);
	}
	auto oneGrid = cpu, observed = cpu;
	const double initial = cpu.Energy(dt) + cpu.ResistorLoss() - cpu.SourceWork();
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(cpu, axis, std::vector<double>{.0047, .0083});
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend*> pointers;
	unsigned int portWorkers = 0, sourceOnlyWorkers = 0;
	for (std::size_t k = 0; k < plan->Domains().size(); ++k)
	{
		const auto &domain = plan->Domains()[k];
		portWorkers += !domain.fields.rootConductance.empty();
		sourceOnlyWorkers += domain.fields.rootConductance.empty() && !domain.fields.sourceEdge.empty();
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(mixed && k % 2 ? "reference" : target.backend));
		Require(bool(backends.back()), "requested backend unavailable");
		(mixed && k % 2 ? SubgridDeviceTarget{} : target).Verify(*backends.back());
		pointers.push_back(backends.back().get());
	}
	Require(portWorkers == 1 && (!soft || sourceOnlyWorkers > 0), "port/source ownership control missing");
	auto domains = CreateCartesianSubgridDomains(plan, pointers);
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> single(GPU_Backend::New(target.backend));
	Require(bool(single), "one-device control unavailable");
	target.Verify(*single);
	auto reference = single->CreateCartesianSubgrid(oneGrid);
	domains->RequireFixedTimestep(dt);
	domains->ReadFields(observed);
	const double scale = std::max({cpu.Energy(dt), cpu.ResistorLoss(), std::abs(cpu.SourceWork())});
	Require(std::abs(observed.ResistorLoss() - cpu.ResistorLoss()) / scale < 2e-6 &&
		std::abs(observed.SourceWork() - cpu.SourceWork()) / scale < 2e-6,
		"initial ledger was duplicated or discarded");
	Reject([&] { domains->Step(dt * .5, 0); });
	for (unsigned int step = 0; step < steps; ++step)
	{
		const double x = (double(step) - 40) / 12;
		const double amplitude = std::exp(-x * x) * (soft && !dc ? x : 1);
		domains->UpdateElectric(dt, amplitude);
		if (step == 0)
		{
			Reject([&] { domains->ReadLedger(); });
		}
		domains->UpdateMagnetic(dt);
		reference->Step(dt, amplitude);
		cpu.Step(dt, amplitude);
		if ((step + 1) % 128 == 0)
		{
			domains->ReadFields(observed);
		}
	}
	domains->ReadFields(observed);
	const auto ledger = domains->ReadLedger();
	Require(ledger[0] == observed.ResistorLoss() && ledger[1] == observed.SourceWork(),
		"explicit ledger read differs from the complete-step observation");
	reference->ReadFields(oneGrid);
	const double denominator = std::max({std::abs(initial), std::abs(cpu.SourceWork()), cpu.ResistorLoss()});
	const double singleError = FieldError(oneGrid, observed), doubleError = FieldError(cpu, observed);
	const double lossError = std::abs(cpu.ResistorLoss() - observed.ResistorLoss()) / denominator;
	const double workError = std::abs(cpu.SourceWork() - observed.SourceWork()) / denominator;
	const double balance = std::abs(observed.Energy(dt) + observed.ResistorLoss() - observed.SourceWork() - initial) / denominator;
	std::cout << "PORT_DOMAINS " << target.backend << " axis=" << axis << " drive=" << drive << " soft=" << soft
		<< " mixed=" << mixed << " dc=" << dc << " steps=" << steps << " single_error=" << singleError
		<< " double_error=" << doubleError << " loss_error=" << lossError << " work_error=" << workError
		<< " balance=" << balance << std::endl;
	Require(singleError < 2e-6 && doubleError < 2e-4, "port-domain field gate failed");
	Require(lossError < 2e-4 && workError < 2e-4 && balance < 2e-4, "port-domain ledger/passivity gate failed");
	Require(observed.ResistorLoss() > 0, "passive ports did not dissipate");
	const auto owner = std::find_if(plan->Domains().begin(), plan->Domains().end(),
		[](const CartesianSubgridDomainData &d) { return !d.fields.rootConductance.empty(); });
	auto invalid = *owner;
	invalid.fields.responseEdge.front() = static_cast<std::uint32_t>(invalid.ownedElectric);
	Reject([&] { CreateCartesianSubgridDomainExecution(*single, invalid); });
	invalid = *owner;
	invalid.fields.rootConductance.front() = 0;
	Reject([&] { CreateCartesianSubgridDomainExecution(*single, invalid); });
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "usage: test_cartesian_subgrid_device_domain_ports [--metal|--cuda|--cuda-two|--cuda-two-reverse]");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		Control(target, 0, false, false, false, false, 4000);
		Control(target, 0, true, false, false, false, 4000);
		Control(target, 1, false, true, false, false, 4000);
		Control(target, 2, true, true, false, false, 4000);
		Control(target, 0, false, true, false, true, 10000);
		if (target.backend != "reference")
		{
			Control(target, 2, true, true, true, false, 4000);
		}
		std::cout << "PASS: logical-domain coupled ports, remote soft work, ledger attribution and fixed-step guards" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

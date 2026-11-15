/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_pml.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
static void Require(bool valid, const char *message)
{
	if (!valid)
		throw std::runtime_error(message);
}

static void Packet(unsigned int normal, int direction, bool refined, bool graded, bool absorb,
				   unsigned int minimum)
{
	const unsigned int en = (normal + 1) % 3, hn = (normal + 2) % 3;
	const double speed = 1 / std::sqrt(epsilon * mu), impedance = std::sqrt(mu / epsilon);
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
		for (unsigned int i = 0; i <= (n == normal ? 80u : 4u); ++i)
			csx.GetGrid()->AddDiscLine(n, i - (graded && n == normal && i % 2 ? .15 : 0));
	std::vector<double> box{1, 1, 1, 3, 3, 3};
	box[normal] = 24;
	box[normal + 3] = 40;
	Grid::Index ratio{1, 1, 1};
	ratio[normal] = refined ? 2 : 1;
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(box, ratio));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "PML geometry failed");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	bc[2 * hn] = bc[2 * hn + 1] = 1;
	std::array<unsigned int, 6> cells{};
	if (absorb)
	{
		bc[2 * normal] = bc[2 * normal + 1] = 3;
		cells[2 * normal] = cells[2 * normal + 1] = 8;
	}
	op->SetBoundaryCondition(bc);
	if (absorb)
		op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), cells));
	Require(op->CalcECOperator() == 0, "PML coefficients failed");
	const double dt = .8 * op->GetTimestep();
	// Rebuild with the same explicit timestep used by the PML recurrence.
	op->SetTimestep(dt);
	Require(op->CalcECOperator() == 0, "PML timestep rebuild failed");
	const double centre = direction > 0 ? .020 : .060, width = .004;
	auto pulse = [&](double p) { return std::exp(-std::pow((p - centre) / width, 2)); };
	op->Grid().SetElectric([&](unsigned int n, const Grid::Point &p)
						   { return n == en ? pulse(p[normal]) : 0.; });
	op->Grid().SetMagnetic(
		[&](unsigned int n, const Grid::Point &p)
		{ return n == hn ? direction * pulse(p[normal] + direction * speed * dt / 2) / impedance : 0.; });
	std::unique_ptr<Engine> eng(op->CreateEngine());
	const double initial = op->Grid().Energy(dt);
	double ledger = 0, peakEnergy = initial, returned = 0;
	Grid::Point probe{.002, .002, .002};
	probe[normal] = direction > 0 ? .016 : .064;
	const unsigned int steps = std::max(minimum, unsigned(5 * .080 / (speed * dt)));
	for (unsigned int step = 0; step < steps; ++step)
	{
		eng->IterateTS(1);
		double work = 0;
		for (std::size_t i = 0; i < eng->GetExtensionCount(); ++i)
			work += dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(eng->GetExtension(i))->BoundaryWorkOut();
		const double energy = op->Grid().Energy(dt);
		ledger = std::max(ledger, std::abs(energy + work - initial) / initial);
		peakEnergy = std::max(peakEnergy, energy);
		if ((step + 1) * dt > 1.3 * .080 / speed)
			returned = std::max(returned, std::abs(op->Grid().Electric(en, probe)));
	}
	const double residual = op->Grid().Energy(dt) / initial;
	std::cout << "PML packet axis " << normal << " sign " << direction << " local " << refined << " graded "
			  << graded << " absorbing " << absorb << " steps " << steps << " return " << returned
			  << " residual " << residual << " ledger " << ledger << " peak-energy " << peakEnergy / initial
			  << std::endl;
	Require(ledger < 1e-9, "PML correction work fails the leapfrog energy identity");
	Require(std::isfinite(residual) && peakEnergy / initial < 1.05, "PML packet unbounded");
	if (absorb)
		Require(returned < .03 && residual < .005, "PML packet reflection/decay failed");
	else
		Require(returned > .5 && residual > .99, "PEC negative control did not retain/reflect the packet");
}

static void Corners(unsigned int minimum, bool multiple = false)
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
		for (unsigned int i = 0; i <= 24; ++i)
			csx.GetGrid()->AddDiscLine(n, i);
	Excitation signal;
	signal.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(
		multiple ? Operator_CartesianSubgrid::NewRegions({{8, 10, 10, 10, 14, 14}, {14, 10, 10, 16, 14, 14}},
														 {{2, 1, 2}, {1, 2, 2}})
				 : Operator_CartesianSubgrid::New({10, 10, 10, 14, 14, 14}));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&signal);
	Require(op->SetGeometryCSX(&csx), "corner geometry failed");
	int bc[6] = {3, 3, 3, 3, 3, 3};
	op->SetBoundaryCondition(bc);
	op->AddExtension(new Operator_Ext_CartesianSubgrid_PML(op.get(), {4, 4, 4, 4, 4, 4}));
	Require(op->CalcECOperator() == 0, "corner coefficients failed");
	const double dt = op->GetTimestep();
	op->Grid().SetElectric(
		[](unsigned int n, const Grid::Point &p)
		{
			const double x = (p[0] - .012) / .002, y = (p[1] - .012) / .002, z = (p[2] - .012) / .002;
			return (n == 1 ? -z : n == 2 ? y : 0) * std::exp(-x * x - y * y - z * z);
		});
	std::unique_ptr<Engine> eng(op->CreateEngine());
	const double initial = op->Grid().Energy(dt);
	double ledger = 0, peak = initial;
	const unsigned int steps = std::max(minimum, 1000u);
	for (unsigned int s = 0; s < steps; ++s)
	{
		eng->IterateTS(1);
		const double energy = op->Grid().Energy(dt);
		const auto *extension = dynamic_cast<Engine_Ext_CartesianSubgrid_PML *>(eng->GetExtension(0));
		ledger = std::max(ledger, std::abs(energy + extension->BoundaryWorkOut() - initial) / initial);
		peak = std::max(peak, energy);
	}
	const double residual = op->Grid().Energy(dt) / initial;
	std::cout << "PML six-face corner multiple " << multiple << " steps " << steps << " residual " << residual
			  << " ledger " << ledger << " peak-energy " << peak / initial << std::endl;
	Require(ledger < 1e-9 && std::isfinite(residual) && peak / initial < 1.1 && residual < .05,
			"PML corner absorption/boundedness/energy ledger failed");
}

int main(int argc, char **argv)
{
	try
	{
		for (unsigned int n = 0; n < 3; ++n)
			for (int sign : {-1, 1})
				Packet(n, sign, false, false, true, 0);
		Packet(0, 1, true, false, true, argc > 1 && std::string(argv[1]) == "--long" ? 100000 : 0);
		Packet(0, -1, true, true, true, 0);
		Packet(0, 1, false, false, false, 0);
		Corners(argc > 1 && std::string(argv[1]) == "--long-corners" ? 10000 : 0);
		Corners(0, true);
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/extensions/engine_ext_cartesiansubgrid_mur.h"
#include "FDTD/extensions/operator_ext_cartesiansubgrid_mur.h"
#include "FDTD/extensions/operator_ext_mur_abc.h"
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
	{
		throw std::runtime_error(message);
	}
}
static void HalfStepContract()
{
	Grid a({6, 6, 6}, {.001, .001, .001}, {{2, 2, 2}, {4, 4, 4}}, epsilon, mu);
	Grid b({6, 6, 6}, {.001, .001, .001}, {{2, 2, 2}, {4, 4, 4}}, epsilon, mu);
	auto seed = [](unsigned int n, const Grid::Point &p)
	{
		return std::sin(1001 * p[0] + 731 * p[1] + 533 * p[2] + n);
	};
	a.SetElectric(seed);
	b.SetElectric(seed);
	const double dt = .7 * a.StableTimestep();
	for (unsigned int i = 0; i < 200; ++i)
	{
		a.Step(dt);
		b.UpdateElectric(dt);
		b.UpdateMagnetic(dt);
	}
	Require(std::abs(a.Energy(dt) / b.Energy(dt) - 1) < 1e-14,
			"split sweeps changed the old full-step update");
	unsigned int rejected = 0;
	try
	{
		b.UpdateMagnetic(dt);
	}
	catch (const std::invalid_argument &)
	{
		++rejected;
	}
	b.UpdateElectric(dt);
	try
	{
		b.UpdateElectric(dt);
	}
	catch (const std::invalid_argument &)
	{
		++rejected;
	}
	try
	{
		b.UpdateMagnetic(dt * .9);
	}
	catch (const std::invalid_argument &)
	{
		++rejected;
	}
	b.UpdateMagnetic(dt);
	Require(rejected == 3, "invalid E/H half-step ordering was accepted");
}
static void Packet(unsigned int normal, int direction, bool refined, unsigned int minimumSteps)
{
	const unsigned int electric = (normal + 1) % 3, magnetic = (normal + 2) % 3;
	const double speed = 1 / std::sqrt(epsilon * mu), impedance = std::sqrt(mu / epsilon);
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= (n == normal ? 80u : 4u); ++i)
		{
			csx.GetGrid()->AddDiscLine(n, i);
		}
	}
	std::vector<double> box = {1, 1, 1, 3, 3, 3};
	box[normal] = 24;
	box[normal + 3] = 40;
	Grid::Index ratio{1, 1, 1};
	ratio[normal] = refined ? 2 : 1;
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New(box, ratio));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(op->SetGeometryCSX(&csx), "Mur geometry setup failed");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	bc[2 * normal] = bc[2 * normal + 1] = 2;
	bc[2 * magnetic] = bc[2 * magnetic + 1] = 1;
	op->SetBoundaryCondition(bc);
	for (unsigned int face = 2 * normal; face <= 2 * normal + 1; ++face)
	{
		op->AddExtension(new Operator_Ext_CartesianSubgrid_Mur(op.get(), face));
	}
	Require(op->CalcECOperator() == 0, "Mur operator preparation failed");
	const double dt = op->GetTimestep();
	const double centre = direction > 0 ? .020 : .060, width = .004;
	auto pulse = [&](double x)
	{
		return std::exp(-std::pow((x - centre) / width, 2));
	};
	auto efield = [&](unsigned int n, const Grid::Point &p)
	{
		return n == electric ? pulse(p[normal]) : 0.;
	};
	auto hfield = [&](unsigned int n, const Grid::Point &p)
	{
		return n == magnetic ? direction * pulse(p[normal] + direction * speed * dt / 2) / impedance : 0.;
	};
	op->Grid().SetElectric(efield);
	op->Grid().SetMagnetic(hfield);
	std::unique_ptr<Engine> eng(op->CreateEngine());
	std::unique_ptr<Operator> yee;
	std::unique_ptr<Engine> yeeEngine;
	if (!refined)
	{
		yee.reset(Operator::New(1));
		yee->SetEngineType(Engine::BASIC);
		yee->SetExcitationSignal(&excitation);
		Require(yee->SetGeometryCSX(&csx), "Yee Mur setup failed");
		yee->SetBoundaryCondition(bc);
		yee->SetTimestep(dt);
		for (unsigned int face = 2 * normal; face <= 2 * normal + 1; ++face)
		{
			auto *ext = new Operator_Ext_Mur_ABC(yee.get());
			ext->SetDirection(normal, face % 2);
			yee->AddExtension(ext);
		}
		Require(yee->CalcECOperator() == 0, "Yee Mur coefficients failed");
		yeeEngine.reset(yee->CreateEngine());
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int x = 0; x < yee->GetNumberOfLines(0); ++x)
			{
				for (unsigned int y = 0; y < yee->GetNumberOfLines(1); ++y)
				{
					for (unsigned int z = 0; z < yee->GetNumberOfLines(2); ++z)
					{
						unsigned int pos[3] = {x, y, z};
						double coord[3];
						if (pos[n] + 1 < yee->GetNumberOfLines(n))
						{
							yee->GetYeeCoords(n, pos, coord, false);
							yeeEngine->SetVolt(
								n, pos,
								efield(n, {coord[0] * .001, coord[1] * .001, coord[2] * .001}) *
									yee->GetEdgeLength(n, pos));
						}
						if (pos[(n + 1) % 3] + 1 < yee->GetNumberOfLines((n + 1) % 3) &&
							pos[(n + 2) % 3] + 1 < yee->GetNumberOfLines((n + 2) % 3))
						{
							yee->GetYeeCoords(n, pos, coord, true);
							yeeEngine->SetCurr(
								n, pos,
								hfield(n, {coord[0] * .001, coord[1] * .001, coord[2] * .001}) *
									yee->GetEdgeLength(n, pos, true));
						}
					}
				}
			}
		}
	}
	const double initial = op->Grid().Energy(dt);
	double residual = 0, returnPeak = 0, yeeError = 0, peakEnergy = initial;
	unsigned int probeIndex[3] = {2, 2, 2};
	probeIndex[normal] = direction > 0 ? 16 : 64;
	const auto probe = op->Position(electric, probeIndex, false);
	const unsigned int steps = std::max(minimumSteps, unsigned(4 * .080 / (speed * dt)));
	for (unsigned int step = 0; step < steps; ++step)
	{
		eng->IterateTS(1);
		double work = 0;
		for (std::size_t i = 0; i < eng->GetExtensionCount(); ++i)
		{
			work += dynamic_cast<Engine_Ext_CartesianSubgrid_Mur *>(eng->GetExtension(i))->BoundaryWorkOut();
		}
		const double energy = op->Grid().Energy(dt);
		residual = std::max(residual, std::abs(energy + work - initial) / initial);
		peakEnergy = std::max(peakEnergy, energy);
		const double field = op->Grid().Electric(electric, probe);
		if ((step + 1) * dt > 1.4 * .080 / speed)
		{
			returnPeak = std::max(returnPeak, std::abs(field));
		}
		if (yeeEngine)
		{
			yeeEngine->IterateTS(1);
			const double reference =
				yeeEngine->GetVolt(electric, probeIndex) / yee->GetEdgeLength(electric, probeIndex);
			yeeError = std::max(yeeError, std::abs(field - reference));
		}
	}
	const double remaining = op->Grid().Energy(dt) / initial;
	std::cout << "Mur packet axis " << normal << " direction " << direction << " refined " << refined
			  << " steps " << steps << " reflected peak " << returnPeak << " remaining energy " << remaining
			  << " energy ledger " << residual << " peak energy ratio " << peakEnergy / initial
			  << " Yee error " << yeeError << '\n';
	Require(residual < 1e-10, "Mur boundary correction work does not balance the leapfrog invariant");
	Require(returnPeak < .03 && remaining < 1e-3 && peakEnergy / initial < 1.01,
			"Mur packet reflection/decay/boundedness failed");
	Require(yeeError < 2e-4, "identity Mur field differs from ordinary Yee");
}
int main(int argc, char **argv)
{
	try
	{
		HalfStepContract();
		for (unsigned int normal = 0; normal < 3; ++normal)
		{
			for (int direction : {-1, 1})
			{
				Packet(normal, direction, false, 0);
			}
		}
		Packet(0, 1, true, argc > 1 && std::string(argv[1]) == "--long" ? 100000 : 0);
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "CSRectGrid.h"
#include "ContinuousStructure.h"
#include "FDTD/cartesian_multigrid_coupling.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesianmultigrid.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

struct Fixture
{
	ContinuousStructure geometry;
	Excitation excitation;
	std::unique_ptr<Operator> op;
	std::unique_ptr<Engine> engine;

	Fixture(Engine::EngineType type, bool subgrid, unsigned int threads = 1, double forcedStep = 0,
			unsigned int refinement = 1)
	{
		CSRectGrid *grid = geometry.GetGrid();
		grid->SetDeltaUnit(1e-3);
		for (int axis = 0; axis < 3; ++axis)
		{
			const int baseCells = axis == 2 ? 8 : 32;
			const int cells = baseCells * refinement;
			for (int i = 0; i <= cells; ++i)
			{
				grid->AddDiscLine(axis, i / double(refinement) - baseCells / 2.);
			}
		}
		excitation.SetupDiracPulse(50e9);
		op.reset(subgrid ? static_cast<Operator *>(Operator_CartesianMultiGrid::New(0, threads))
						 : Operator::New(threads));
		op->SetEngineType(type);
		op->SetExcitationSignal(&excitation);
		Require(op->SetGeometryCSX(&geometry), "geometry setup failed");
		Require(grid->GetQtyLines(0) == 32 * refinement + 1 && grid->GetQtyLines(1) == 32 * refinement + 1,
				"shared grid was not restored");
		int boundaries[6] = {0, 0, 0, 0, 0, 0};
		op->SetBoundaryCondition(boundaries);
		op->SetTimeStepMethod(1);
		op->SetTimestepFactor(forcedStep > 0 ? 1 : .8);
		if (forcedStep > 0)
		{
			op->SetTimestep(forcedStep);
		}
		Require(op->CalcECOperator() == 0, "operator build failed");
		if (subgrid)
		{
			auto *multi = static_cast<Operator_CartesianMultiGrid *>(op.get());
			Require(multi->GetCoarseOperator()->GetTimestep() == op->GetTimestep(), "unequal timesteps");
			Require(multi->GetCoupledStepNormBound() < 2, "coupled-curl timestep bound failed");
			const auto &ip = multi->GetInterpolation();
			for (int direction = 0; direction < 2; ++direction)
			{
				for (unsigned int y = 0; y < 32 * refinement + 1; ++y)
				{
					Require(ip.pos_v_2p[direction][y] < 16 * refinement + 1 &&
								ip.pos_v_2pp[direction][y] < 16 * refinement + 1,
							"voltage interpolation outside child");
					Require(ip.pos_i_2p[direction][y] < 16 * refinement + 1 &&
								ip.pos_i_2pp[direction][y] < 16 * refinement + 1,
							"current interpolation outside child");
				}
			}
		}
		engine.reset(op->CreateEngine());
		Seed(engine.get(), op.get());
		if (subgrid)
		{
			Seed(engine->GetSubGridEngine(), static_cast<Operator_CartesianMultiGrid *>(op.get())->GetCoarseOperator());
		}
	}

	static void Seed(Engine *engine, Operator *op)
	{
		const double pi = std::acos(-1.);
		for (unsigned int x = 0; x < op->GetNumberOfLines(0); ++x)
		{
			for (unsigned int y = 0; y < op->GetNumberOfLines(1); ++y)
			{
				const double field =
					std::sin(pi * (op->GetDiscLine(0, x) + 16) / 32) * std::sin(pi * (op->GetDiscLine(1, y) + 16) / 32);
				for (unsigned int z = 0; z + 1 < op->GetNumberOfLines(2); ++z)
				{
					engine->SetVolt(2, x, y, z, field * op->GetDiscDelta(2, z) * op->GetGridDelta());
				}
			}
		}
	}

	double Energy() const
	{
		double sum = 0;
		for (unsigned int x = 0; x < 33; ++x)
		{
			for (unsigned int y = 0; y < 33; ++y)
			{
				for (unsigned int z = 0; z < 9; ++z)
				{
					for (unsigned int n = 0; n < 3; ++n)
					{
						const unsigned int pos[3] = {x, y, z};
						const double e = engine->GetVolt(n, x, y, z) / op->GetEdgeLength(n, pos);
						const double h = engine->GetCurr(n, x, y, z) / op->GetEdgeLength(n, pos, true);
						Require(std::isfinite(e) && std::isfinite(h), "nonfinite cavity field");
						sum += e * e + 376.730313668 * 376.730313668 * h * h;
					}
				}
			}
		}
		return sum;
	}
};

static void Reject(const std::function<void()> &action)
{
	bool rejected = false;
	try
	{
		action();
	}
	catch (const std::exception &)
	{
		rejected = true;
	}
	Require(rejected, "invalid subgrid was silently accepted");
}

static void CompareEngines(Engine::EngineType type, unsigned int threads)
{
	Fixture reference(Engine::BASIC, true);
	Fixture candidate(type, true, threads);
	reference.engine->IterateTS(256);
	candidate.engine->IterateTS(256);
	double error = 0, peak = 0, currentError = 0, currentPeak = 0;
	for (unsigned int x = 0; x < 33; ++x)
	{
		for (unsigned int y = 0; y < 33; ++y)
		{
			for (unsigned int z = 0; z < 9; ++z)
			{
				for (unsigned int n = 0; n < 3; ++n)
				{
					const double a = reference.engine->GetVolt(n, x, y, z);
					const double b = candidate.engine->GetVolt(n, x, y, z);
					peak = std::max(peak, std::abs(a));
					error = std::max(error, std::abs(a - b));
					const double ha = reference.engine->GetCurr(n, x, y, z);
					const double hb = candidate.engine->GetCurr(n, x, y, z);
					currentPeak = std::max(currentPeak, std::abs(ha));
					currentError = std::max(currentError, std::abs(ha - hb));
				}
			}
		}
	}
	std::cout << "engine " << type << ", threads " << threads << ": E relative error " << error / peak << '\n';
	Require(peak > 0 && error / peak < 2e-4, "engine E comparison failed");
	std::cout << "engine H relative error " << currentError / currentPeak << '\n';
	Require(currentPeak > 0 && currentError / currentPeak < 2e-4, "engine H comparison failed");
}

static void CheckPowerAdjoint()
{
	Fixture fixture(Engine::BASIC, true);
	auto *op = static_cast<Operator_CartesianMultiGrid *>(fixture.op.get());
	Engine *fine = fixture.engine.get();
	Engine *coarse = fine->GetSubGridEngine();
	const unsigned int lines[3] = {33, 33, 9};
	CartesianMultiGridCoupling coupling(fine, coarse, op->GetInterpolation(), lines);
	std::mt19937 generator(2026);
	std::uniform_real_distribution<float> value(-1, 1);
	for (unsigned int y = 0; y < 32; ++y)
	{
		for (unsigned int z = 0; z < 8; ++z)
		{
			fine->SetVolt(1, 16, y, z, value(generator));
			fine->SetVolt(2, 16, y, z, y == 0 ? 0 : value(generator));
		}
	}
	for (unsigned int y = 0; y < 16; ++y)
	{
		for (unsigned int z = 0; z < 8; ++z)
		{
			coarse->SetCurr(1, 15, y, z, value(generator));
			coarse->SetCurr(2, 15, y, z, value(generator));
		}
	}
	coupling.SyncVoltages();
	coupling.SyncCurrents();
	double a = 0, b = 0, scale = 0;
	for (unsigned int y = 0; y < 32; ++y)
	{
		for (unsigned int z = 0; z < 8; ++z)
		{
			const double term = double(fine->GetVolt(1, 16, y, z)) * fine->GetCurr(2, 15, y, z) -
								double(fine->GetVolt(2, 16, y, z)) * fine->GetCurr(1, 15, y, z);
			a += term;
			scale += std::abs(term);
		}
	}
	for (unsigned int y = 0; y < 16; ++y)
	{
		for (unsigned int z = 0; z < 8; ++z)
		{
			b += double(coarse->GetVolt(1, 16, y, z)) * coarse->GetCurr(2, 15, y, z) -
				 double(coarse->GetVolt(2, 16, y, z)) * coarse->GetCurr(1, 15, y, z);
		}
	}
	std::cout << "interface power mismatch / absolute sum " << std::abs(a - b) / scale << '\n';
	Require(scale > 0 && std::abs(a - b) < 1e-6 * scale, "interface transfer is not power-adjoint");
}

static double PacketError(unsigned int refinement, unsigned int component)
{
	Fixture candidate(Engine::BASIC, true, 1, 0, refinement);
	Fixture reference(Engine::BASIC, false, 1, candidate.op->GetTimestep(), refinement);
	auto seed = [component](Engine *engine, Operator *op)
	{
		for (unsigned int x = 0; x < op->GetNumberOfLines(0); ++x)
		{
			for (unsigned int y = 0; y < op->GetNumberOfLines(1); ++y)
			{
				const double field =
					std::exp(-std::pow((op->GetDiscLine(0, x) + 8) / 2., 2)) *
					(component == 2 ? std::sin(std::acos(-1.) * (op->GetDiscLine(1, y) + 16) / 32.) : 1.);
				for (unsigned int z = 0; z < op->GetNumberOfLines(2); ++z)
				{
					for (unsigned int n = 0; n < 3; ++n)
					{
						engine->SetVolt(
							n, x, y, z,
							n == component &&
									(component == 2 ? z + 1 < op->GetNumberOfLines(2) : y + 1 < op->GetNumberOfLines(1))
								? field *
									  (component == 2 ? 1.
													  : std::sin(std::acos(-1.) * (op->GetDiscLine(2, z) + 4) / 8.)) *
									  op->GetDiscDelta(component, component == 2 ? z : y) * op->GetGridDelta()
								: 0);
						engine->SetCurr(n, x, y, z, 0);
					}
				}
			}
		}
	};
	seed(candidate.engine.get(), candidate.op.get());
	seed(candidate.engine->GetSubGridEngine(),
		 static_cast<Operator_CartesianMultiGrid *>(candidate.op.get())->GetCoarseOperator());
	seed(reference.engine.get(), reference.op.get());
	double error = 0, peak = 0;
	for (unsigned int step = 0; step < 60 * refinement; ++step)
	{
		candidate.engine->IterateTS(1);
		reference.engine->IterateTS(1);
		const double a = reference.engine->GetVolt(component, 24 * refinement, 16 * refinement, 4 * refinement);
		const double b = candidate.engine->GetVolt(component, 24 * refinement, 16 * refinement, 4 * refinement);
		peak = std::max(peak, std::abs(a));
		error = std::max(error, std::abs(a - b));
	}
	std::cout << "packet refinement " << refinement << ", E component " << component << ": peak relative error "
			  << error / peak << '\n';
	Require(peak > 1e-6, "no transmitted pulse in reference");
	return error / peak;
}

static void CheckPacketTransmission()
{
	for (unsigned int component : {1u, 2u})
	{
		const double coarse = PacketError(1, component);
		const double refined = PacketError(2, component);
		const double finest = PacketError(4, component);
		Require((refined < coarse && finest < refined) || coarse < 1e-6,
				"interface packet error does not decrease on refinement");
		Require(finest < .02, "resolved transmitted pulse differs by more than two percent");
	}
}

static void CheckRandomStability()
{
	Fixture fixture(Engine::BASIC, true);
	std::mt19937 generator(42);
	std::uniform_real_distribution<float> random(-1, 1);
	auto seed = [&](Engine *engine, Operator *op)
	{
		for (unsigned int x = 1; x + 1 < op->GetNumberOfLines(0); ++x)
		{
			for (unsigned int y = 1; y + 1 < op->GetNumberOfLines(1); ++y)
			{
				for (unsigned int z = 1; z + 1 < op->GetNumberOfLines(2); ++z)
				{
					const unsigned int pos[3] = {x, y, z};
					for (unsigned int n = 0; n < 3; ++n)
					{
						engine->SetVolt(n, x, y, z, random(generator) * op->GetEdgeLength(n, pos));
						engine->SetCurr(n, x, y, z,
										random(generator) * op->GetEdgeLength(n, pos, true) / 376.730313668);
					}
				}
			}
		}
	};
	seed(fixture.engine.get(), fixture.op.get());
	seed(fixture.engine->GetSubGridEngine(),
		 static_cast<Operator_CartesianMultiGrid *>(fixture.op.get())->GetCoarseOperator());
	const double initial = fixture.Energy();
	double maximum = initial;
	for (unsigned int step = 0; step < 20000; step += 100)
	{
		fixture.engine->IterateTS(100);
		maximum = std::max(maximum, fixture.Energy());
		Require(maximum < 5 * initial, "random three-component fields are unstable");
	}
	std::cout << "random E/H cavity: max field norm / initial " << maximum / initial << '\n';
}

int main(int argc, char **argv)
{
	try
	{
		Reject([] { std::unique_ptr<Operator_CartesianMultiGrid> op(Operator_CartesianMultiGrid::New(NAN)); });
		Reject(
			[]
			{
				std::unique_ptr<Operator_CartesianMultiGrid> op(Operator_CartesianMultiGrid::New(0));
				op->SetEngineType(Engine::GPU);
			});
		Reject(
			[]
			{
				std::unique_ptr<Operator_CartesianMultiGrid> op(Operator_CartesianMultiGrid::New(0));
				int boundaries[6] = {0, 0, 0, 0, 0, 3};
				op->SetBoundaryCondition(boundaries);
			});
		CompareEngines(Engine::BASIC, 2);
		CompareEngines(Engine::SSE, 2);
		CompareEngines(Engine::CPU, 2);
		CheckPowerAdjoint();
		CheckPacketTransmission();
		CheckRandomStability();
		Fixture fine(Engine::BASIC, false);
		Fixture refined(Engine::BASIC, true);
		const double initial = refined.Energy();
		double maximum = initial;
		std::vector<double> crossing;
		std::vector<double> fineCrossing;
		double previous = refined.engine->GetVolt(2, 8, 8, 4);
		double finePrevious = fine.engine->GetVolt(2, 8, 8, 4);
		const unsigned int steps = argc == 2 && std::string(argv[1]) == "--long" ? 1000000 : 30000;
		for (unsigned int step = 1; step <= steps; ++step)
		{
			refined.engine->IterateTS(1);
			fine.engine->IterateTS(1);
			const double value = refined.engine->GetVolt(2, 8, 8, 4);
			const double fineValue = fine.engine->GetVolt(2, 8, 8, 4);
			if (previous <= 0 && value > 0)
			{
				crossing.push_back(step - value / (value - previous));
			}
			previous = value;
			if (finePrevious <= 0 && fineValue > 0)
			{
				fineCrossing.push_back(step - fineValue / (fineValue - finePrevious));
			}
			finePrevious = fineValue;
			if (step % 100 == 0)
			{
				maximum = std::max(maximum, refined.Energy());
				if (maximum >= 3 * initial)
				{
					std::cerr << "unstable step " << step << ", energy/initial " << maximum / initial << '\n';
				}
				Require(maximum < 3 * initial, "late-time cavity instability");
			}
		}
		Require(crossing.size() > 20, "no transmitted cavity oscillation");
		Require(fineCrossing.size() > 20, "no fine reference cavity oscillation");
		const double frequency =
			(crossing.size() - 1) / ((crossing.back() - crossing.front()) * refined.op->GetTimestep());
		const double exact = 299792458. / 2 * std::sqrt(2.) / .032;
		const double fineFrequency =
			(fineCrossing.size() - 1) / ((fineCrossing.back() - fineCrossing.front()) * fine.op->GetTimestep());
		std::cout << steps << "-step cavity: energy max/initial " << maximum / initial << ", frequency " << frequency
				  << ", analytic " << exact << '\n';
		std::cout << "uniform-fine frequency " << fineFrequency << ", relative difference "
				  << frequency / fineFrequency - 1 << '\n';
		Require(std::abs(frequency / exact - 1) < .02, "cavity frequency differs by more than two percent");
		Require(std::abs(frequency / fineFrequency - 1) < .01,
				"cavity frequency differs from uniform fine by more than one percent");
		Require(refined.Energy() > .1 * initial, "coupling erased the cavity energy");
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
	return 0;
}

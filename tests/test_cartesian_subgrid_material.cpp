/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "CSPrimBox.h"
#include "CSPropMaterial.h"
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/operator_cartesiansubgrid.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12;
constexpr double mu = 1.2566370614359173e-6;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

static Grid::Axes Lines(unsigned int cells)
{
	Grid::Axes lines;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= cells; ++i)
		{
			lines[n].push_back(.012 * i / cells);
		}
	}
	return lines;
}

static std::vector<Grid::Region> Regions(unsigned int cells)
{
	return {{{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, {2, 2, 2}}};
}

static double Seed(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

static void Invariant(unsigned int steps)
{
	// Both dielectric and permeability discontinuities cross refinement faces.
	auto medium = [](const Grid::Point &lo, const Grid::Point &hi)
	{
		return Grid::Medium{epsilon * ((lo[0] + hi[0]) / 2 < .006 ? 1 : 3.556),
							mu * ((lo[1] + hi[1]) / 2 < .006 ? 1 : 1.7)};
	};
	Grid grid(Lines(6), epsilon, mu, Regions(6), false, medium);
	grid.SetElectric(Seed);
	grid.SetMagnetic(
		[](unsigned int n, const Grid::Point &p)
		{
			return Seed(n, p) / 376.730313668;
		});
	const double dt = .8 * grid.StableTimestep(), initial = grid.Energy(dt);
	double drift = 0;
	for (unsigned int i = 0; i < steps; ++i)
	{
		grid.Step(dt);
		if (i % 100 == 0)
		{
			drift = std::max(drift, std::abs(grid.Energy(dt) / initial - 1));
		}
	}
	Require(std::isfinite(drift) && drift < 1e-8, "heterogeneous lossless invariant failed");
	Require(grid.CurlDivergenceResidual(Seed) < 1e-13, "material coefficients broke div(curl)");
	std::cout << "heterogeneous invariant " << steps << " steps: " << drift << '\n';
	const Grid::Point boundary{.00317, .006, .00547};
	const double value = grid.Magnetic(1, boundary);
	for (double offset : {-2e-18, 2e-18})
	{
		auto point = boundary;
		point[1] += offset;
		Require(std::abs(grid.Magnetic(1, point) - value) < 1e-10 * std::max(std::abs(value), 1e-30),
				"coordinate roundoff changed the selected material side");
	}
	auto positive = boundary;
	positive[1] += 1e-11;
	Require(std::abs(grid.Magnetic(1, positive) - value) < 1e-6 * std::max(std::abs(value), 1e-30),
			"interface sampling does not select the positive-coordinate leaf");
	for (double value :
		 {0., -epsilon, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
	{
		for (bool electric : {true, false})
		{
			bool rejected = false;
			try
			{
				Grid bad(Lines(6), epsilon, mu, Regions(6), false,
						 [&](const Grid::Point &, const Grid::Point &)
						 {
							 return Grid::Medium{electric ? value : epsilon, electric ? mu : value};
						 });
			}
			catch (const std::invalid_argument &)
			{
				rejected = true;
			}
			Require(rejected, "invalid leaf material accepted");
		}
	}
}

static void ConstantMediumEquivalence()
{
	Grid original(Lines(6), 2.25 * epsilon, 1.3 * mu, Regions(6));
	Grid mapped(Lines(6), epsilon, mu, Regions(6), false,
				[](const Grid::Point &, const Grid::Point &)
				{
					return Grid::Medium{2.25 * epsilon, 1.3 * mu};
				});
	Require(std::abs(original.StableTimestep() / mapped.StableTimestep() - 1) < 1e-13,
			"constant material map changed timestep");
	original.SetElectric(Seed);
	mapped.SetElectric(Seed);
	original.SetMagnetic(Seed);
	mapped.SetMagnetic(Seed);
	const double dt = .8 * original.StableTimestep();
	for (unsigned int i = 0; i < 300; ++i)
	{
		original.Step(dt);
		mapped.Step(dt);
	}
	Require(std::abs(original.Energy(dt) / mapped.Energy(dt) - 1) < 1e-12,
			"constant material map changed field histories");
	for (unsigned int n = 0; n < 3; ++n)
	{
		const Grid::Point p{.00317, .00431, .00547};
		Require(std::abs(original.Electric(n, p) - mapped.Electric(n, p)) < 1e-12 &&
					std::abs(original.Magnetic(n, p) - mapped.Magnetic(n, p)) < 1e-12,
				"constant material map differs from homogeneous controls");
	}
}

// TE cavity split into two dielectrics at x=L/2. E_z and its normal
// derivative are continuous; k1*cot(k1*a)+k2*cot(k2*b)=0, mu identical.
static double LayeredCavity(unsigned int cells)
{
	const double length = .012, pi = std::acos(-1.), transverse = pi / length;
	auto derivative = [&](double wave, double eps)
	{
		const double square = eps * wave * wave - transverse * transverse;
		if (square > 0)
		{
			const double k = std::sqrt(square);
			return k / std::tan(k * length / 2);
		}
		if (square < 0)
		{
			const double k = std::sqrt(-square);
			return k / std::tanh(k * length / 2);
		}
		return 2 / length;
	};
	double lo = transverse / 1.5, hi = std::sqrt(5.) * transverse / 1.5 * (1 - 1e-10);
	for (unsigned int i = 0; i < 100; ++i)
	{
		const double mid = (lo + hi) / 2;
		if (derivative(mid, 1) + derivative(mid, 2.25) > 0)
		{
			lo = mid;
		}
		else
		{
			hi = mid;
		}
	}
	const double wave = (lo + hi) / 2;
	const double k1 = std::sqrt(wave * wave - transverse * transverse);
	const double k2 = std::sqrt(2.25 * wave * wave - transverse * transverse);
	const double amplitude = std::sin(k1 * length / 2) / std::sin(k2 * length / 2);
	const double exact = wave / (2 * pi * std::sqrt(epsilon * mu));
	Grid grid(Lines(cells), epsilon, mu, Regions(cells), false,
			  [&](const Grid::Point &a, const Grid::Point &b)
			  {
				  return Grid::Medium{epsilon * ((a[0] + b[0]) / 2 < length / 2 ? 1 : 2.25), mu};
			  });
	grid.SetElectric(
		[&](unsigned int n, const Grid::Point &p)
		{
			return n == 2 ? std::sin(transverse * p[1]) * (p[0] < length / 2
															   ? std::sin(k1 * p[0])
															   : amplitude * std::sin(k2 * (length - p[0])))
						  : 0.;
		});
	const Grid::Point probe{.00317, .00431, .00547};
	const double dt = .8 * grid.StableTimestep();
	double previous = grid.Electric(2, probe);
	std::vector<double> crossing;
	for (unsigned int step = 1; step <= 4000 * cells / 6; ++step)
	{
		grid.Step(dt);
		const double value = grid.Electric(2, probe);
		if (previous <= 0 && value > 0)
		{
			crossing.push_back(step - value / (value - previous));
		}
		previous = value;
	}
	Require(crossing.size() > 15, "no layered cavity mode observed");
	const double frequency = (crossing.size() - 1) / ((crossing.back() - crossing.front()) * dt);
	const double error = frequency / exact - 1;
	std::cout << "layered TE cavity base " << cells << ": frequency " << frequency << " analytic " << exact
			  << " relative error " << error << '\n';
	return std::abs(error);
}

static void NativeOccupancy()
{
	// Nonzero drawing origin, millimetre units, overlapping priorities and a
	// dielectric crossing the coarse/fine boundary in all three directions.
	ContinuousStructure csx;
	const Grid::Point origin{31, -23, 17};
	csx.GetGrid()->SetDeltaUnit(1e-3);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 6; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, origin[n] + i);
		}
	}
	auto add = [&](double eps, double permeability, unsigned int priority, Grid::Point lo, Grid::Point hi)
	{
		auto *material = new CSPropMaterial(csx.GetParameterSet());
		material->SetEpsilon(eps);
		material->SetMue(permeability);
		csx.AddProperty(material);
		auto *box = new CSPrimBox(csx.GetParameterSet(), material);
		box->SetPriority(priority);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box->SetCoord(2 * n, origin[n] + lo[n]);
			box->SetCoord(2 * n + 1, origin[n] + hi[n]);
		}
	};
	// Insert the higher priority first: insertion order is not the priority.
	add(6, 1.7, 9, {3, 2, 2}, {4, 4, 4});
	add(3.556, 1, 2, {3, 0, 0}, {6, 6, 6});
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(
		Operator_CartesianSubgrid::New({32, -22, 18, 36, -18, 22}, {2, 1, 4}));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(op->SetGeometryCSX(&csx), "native material geometry rejected");
	Require(op->CalcECOperator() == 0, "native material operator failed");
	const auto &grid = op->Grid();
	unsigned int occupied = 0, priority = 0;
	for (std::size_t i = 0; i < grid.CellCount(); ++i)
	{
		const auto lo = grid.LeafStart(i), length = grid.LeafLengths(i);
		Grid::Point centre;
		for (unsigned int n = 0; n < 3; ++n)
		{
			centre[n] = (lo[n] + length[n] / 2) * 1000;
		}
		const bool high = centre[0] > 3 && centre[0] < 4 && centre[1] > 2 && centre[1] < 4 && centre[2] > 2 &&
						  centre[2] < 4;
		const double wantEps = high ? 6 : centre[0] > 3 ? 3.556 : 1;
		const double wantMu = high ? 1.7 : 1;
		const auto value = grid.LeafMedium(i);
		Require(std::abs(value.epsilon / epsilon - wantEps) < 1e-8 &&
					std::abs(value.permeability / mu - wantMu) < 1e-8,
				"native priority, units or material occupancy wrong");
		occupied += wantEps != 1;
		priority += high;
	}
	Require(occupied > priority && priority > 0, "native material test did not exercise both priorities");
	std::cout << "native material occupancy checked " << grid.CellCount() << " leaves, " << occupied
			  << " material and " << priority << " high-priority leaves\n";
}

int main(int argc, char **argv)
{
	try
	{
		ConstantMediumEquivalence();
		const double coarse = LayeredCavity(6), fine = LayeredCavity(12);
		Require(fine < .01 && fine < .8 * coarse, "layered cavity failed frequency convergence");
		NativeOccupancy();
		Invariant(argc == 2 && std::string(argv[1]) == "--long" ? 100000 : 4000);
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	return 0;
}

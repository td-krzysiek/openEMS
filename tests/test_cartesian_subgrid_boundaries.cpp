/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
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
	for (auto &axis : lines)
	{
		for (unsigned int i = 0; i <= cells; ++i)
		{
			axis.push_back(.012 * i / cells);
		}
	}
	return lines;
}
static std::vector<Grid::Region> Regions(unsigned int cells)
{
	return {{{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, {2, 2, 2}}};
}
static Grid::BoundaryConditions PEC()
{
	return {Grid::BoundaryType::PEC, Grid::BoundaryType::PEC, Grid::BoundaryType::PEC,
			Grid::BoundaryType::PEC, Grid::BoundaryType::PEC, Grid::BoundaryType::PEC};
}
static void TracesAndGuards()
{
	const auto lines = Lines(6);
	Grid reference(lines, epsilon, mu, Regions(6));
	for (unsigned int face = 0; face < 6; ++face)
	{
		auto boundary = PEC();
		boundary[face] = Grid::BoundaryType::PMC;
		Grid grid(lines, epsilon, mu, Regions(6), false, {}, {}, boundary);
		grid.SetElectric(
			[](unsigned int, const Grid::Point &)
			{
				return 1.;
			});
		const unsigned int normal = face / 2;
		Grid::Point p{.006, .006, .006};
		p[normal] = face % 2 ? .012 : 0;
		Require(std::abs(grid.Electric((normal + 1) % 3, p) - 1) < 1e-12,
				"PMC incorrectly eliminated a tangential E trace");
		p[normal] = face % 2 ? 0 : .012;
		Require(std::abs(grid.Electric((normal + 1) % 3, p)) < 1e-12, "opposite PEC trace was released");
		Require(grid.EdgeCount() > reference.EdgeCount(), "PMC added no boundary degrees of freedom");
	}
	bool rejected = false;
	try
	{
		auto boundary = PEC();
		boundary[2] = static_cast<Grid::BoundaryType>(9);
		Grid bad(lines, epsilon, mu, Regions(6), false, {}, {}, boundary);
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "unknown boundary condition accepted");
	std::cout << "all six PMC/PEC tangential trace selections and invalid boundary: PASS\n";
}
static void Invariant(unsigned int steps)
{
	auto boundary = PEC();
	for (unsigned int i = 0; i < 6; ++i)
	{
		boundary[i] = Grid::BoundaryType::PMC;
	}
	Grid grid(Lines(8), epsilon, mu,
			  {{{{2, 2, 2}, {4, 4, 4}}, {2, 1, 4}}, {{{5, 2, 2}, {7, 4, 4}}, {1, 2, 1}}}, false, {}, {},
			  boundary);
	auto seed = [](unsigned int n, const Grid::Point &p)
	{
		return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n);
	};
	grid.SetElectric(seed);
	grid.SetMagnetic(
		[&](unsigned int n, const Grid::Point &p)
		{
			return seed(n, p) / 376.730313668;
		});
	const double dt = .8 * grid.StableTimestep(), initial = grid.Energy(dt);
	double drift = 0;
	for (unsigned int step = 0; step < steps; ++step)
	{
		grid.Step(dt);
		if (step % 100 == 0)
		{
			drift = std::max(drift, std::abs(grid.Energy(dt) / initial - 1));
		}
	}
	Require(std::isfinite(drift) && drift < 1e-9, "PMC lossless invariant drift");
	Require(grid.CurlDivergenceResidual(seed) < 1e-13, "PMC broke div(curl)");
	std::cout << "PMC multi-region " << steps << " steps invariant drift " << drift << '\n';
}
static double Cavity(unsigned int cells, unsigned int component, bool allPMC, bool reverse)
{
	const double pi = std::acos(-1.), length = .012;
	const unsigned int a = (component + 1) % 3, b = (component + 2) % 3;
	auto boundaries = PEC();
	if (allPMC)
	{
		for (auto &boundary : boundaries)
		{
			boundary = Grid::BoundaryType::PMC;
		}
	}
	else
	{
		boundaries[2 * a + (reverse ? 0 : 1)] = Grid::BoundaryType::PMC;
		boundaries[2 * b + (reverse ? 1 : 0)] = Grid::BoundaryType::PMC;
	}
	Grid grid(Lines(cells), epsilon, mu, Regions(cells), false, {}, {}, boundaries);
	grid.SetElectric(
		[&](unsigned int n, const Grid::Point &p)
		{
			if (allPMC)
			{
				if (n == a)
				{
					return std::sin(pi * p[a] / length) * std::cos(pi * p[b] / length);
				}
				if (n == b)
				{
					return -std::cos(pi * p[a] / length) * std::sin(pi * p[b] / length);
				}
				return 0.;
			}
			const double x = pi * p[a] / (2 * length), y = pi * p[b] / (2 * length);
			return n == component ? (reverse ? std::cos(x) * std::sin(y) : std::sin(x) * std::cos(y)) : 0.;
		});
	const double dt = .8 * grid.StableTimestep();
	const Grid::Point receiver{.00317, .00431, .00547};
	const unsigned int observed = allPMC ? a : component;
	double previous = grid.Electric(observed, receiver);
	std::vector<double> crossings;
	for (unsigned int step = 1; step <= 3000 * cells / 6; ++step)
	{
		grid.Step(dt);
		const double value = grid.Electric(observed, receiver);
		if (previous <= 0 && value > 0)
		{
			crossings.push_back(step - value / (value - previous));
		}
		previous = value;
	}
	Require(crossings.size() > 10, "no mixed-boundary cavity mode observed");
	const double frequency = (crossings.size() - 1) / ((crossings.back() - crossings.front()) * dt);
	const double exact = std::sqrt(allPMC ? 2. : .5) / (2 * length * std::sqrt(epsilon * mu));
	const double error = std::abs(frequency / exact - 1);
	std::cout << (allPMC ? "all PMC" : "mixed PEC/PMC") << " component " << component << " reverse "
			  << reverse << " cells " << cells << " frequency relative error " << error << '\n';
	return error;
}
int main(int argc, char **argv)
{
	try
	{
		TracesAndGuards();
		Invariant(argc > 1 && std::string(argv[1]) == "--long" ? 100000 : 4000);
		for (bool allPMC : {false, true})
		{
			for (unsigned int component = 0; component < 3; ++component)
			{
				const bool reverse = component == 1;
				const double coarse = Cavity(6, component, allPMC, reverse);
				const double fine = Cavity(12, component, allPMC, reverse);
				Require(fine < coarse * .6 && fine < .01,
						"PMC cavity does not converge toward analytic mode");
			}
		}
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

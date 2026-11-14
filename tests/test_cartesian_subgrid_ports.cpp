/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
constexpr double eps = 8.85418781762e-12, mu = 1.2566370614359173e-6;
static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}
static Grid MakeGrid()
{
	return Grid({6, 6, 6}, {.002, .002, .002}, .0 + eps, mu, {{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}});
}
static void Integrals()
{
	auto grid = MakeGrid();
	grid.SetElectric(
		[](unsigned int n, const Grid::Point &p)
		{
			return 2 + p[(n + 1) % 3] * 100 + p[(n + 2) % 3] * 200;
		});
	const Grid::Point a{.006, .006, .002}, b{.006, .006, .010}, c{.006, .006, .006};
	const auto full = grid.VoltagePath(a, b), first = grid.VoltagePath(a, c), second = grid.VoltagePath(c, b);
	const double value = grid.Evaluate(full);
	Require(std::abs(value - 3.8 * .008) < 1e-14, "voltage is not an exact local-cell integral");
	Require(std::abs(value + grid.Evaluate(grid.VoltagePath(b, a))) < 1e-14, "voltage orientation failed");
	Require(std::abs(value - grid.Evaluate(first) - grid.Evaluate(second)) < 1e-14,
			"voltage subdivision failed");
	grid.SetMagnetic(
		[](unsigned int n, const Grid::Point &p)
		{
			return n == 0 ? -p[1] : n == 1 ? p[0] : 0;
		});
	const auto loop = grid.CurrentContour(2, {.003, .003, .006}, {.009, .009, .006});
	Require(std::abs(grid.Evaluate(loop) - 2 * .006 * .006) < 1e-15,
			"current contour orientation/units failed");
	const auto bounds = grid.DualBounds(0, {.006, .006, .006});
	Require(std::abs(bounds.first - .0055) < 1e-15 && std::abs(bounds.second - .0065) < 1e-15,
			"dual contour used base rather than local spacing");
	const auto rounded = grid.DualBounds(0, {std::nextafter(.006, 0.), .006, .006});
	Require(std::abs(rounded.first - bounds.first) < 1e-15 &&
				std::abs(rounded.second - bounds.second) < 1e-15,
			"coordinate roundoff collapsed a dual contour");
	std::cout << "exact U path across coarse/fine cells and H circulation: PASS\n";
}
static void Power(unsigned int steps)
{
	auto grid = MakeGrid();
	const double dt = .7 * grid.StableTimestep();
	grid.SetResistors({{{.002, .006, .004}, {.002, .006, .006}, 50, .02},
					   {{.006, .006, .005}, {.006, .006, .006}, 50, 0},
					   {{.004, .006, .002}, {.004, .006, .010}, 75, 0}},
					  dt);
	double peak = 0, error = 0, passiveIncrease = 0, previous = 0;
	for (unsigned int step = 0; step < steps; ++step)
	{
		const double drive =
			step < 400 ? std::sin(.07 * step) * std::pow(std::sin(std::acos(-1.) * step / 400), 2) : 0;
		grid.Step(dt, drive);
		const double energy = grid.Energy(dt);
		peak = std::max(peak, grid.SourceWork());
		error = std::max(error, std::abs(energy + grid.ResistorLoss() - grid.SourceWork()));
		if (step > 400)
		{
			passiveIncrease = std::max(passiveIncrease, energy - previous);
		}
		previous = energy;
	}
	Require(peak > 0 && error / peak < 1e-10 && passiveIncrease / peak < 1e-12,
			"resistor passivity or injected/dissipated/stored energy balance failed");
	Require(grid.ResistorLoss() > 0, "resistors did not dissipate power");
	std::cout << steps << " steps, normalized energy residual " << error / peak << ", passive increase "
			  << passiveIncrease / peak << '\n';
}
static void SoftSourcePower()
{
	auto grid = MakeGrid();
	const double dt = .7 * grid.StableTimestep();
	grid.SetResistors({{{.006, .006, .005}, {.006, .006, .006}, 50, 0}}, dt);
	grid.SetSources({{grid.EdgeCount() / 2, 1}});
	double peak = 0, error = 0;
	for (unsigned int step = 0; step < 1000; ++step)
	{
		grid.Step(dt, step < 200 ? std::sin(.23 * step) : 0);
		peak = std::max(peak, grid.SourceWork());
		error = std::max(error, std::abs(grid.Energy(dt) + grid.ResistorLoss() - grid.SourceWork()));
	}
	Require(peak > 0 && error / peak < 1e-10, "soft source plus resistor lost energy accounting");
	std::cout << "soft E source + resistor energy residual " << error / peak << '\n';
}
static void ParallelAndGuards()
{
	auto one = MakeGrid(), two = MakeGrid();
	const double dt = .6 * one.StableTimestep();
	const Grid::Point a{.006, .006, .005}, b{.006, .006, .006};
	one.SetResistors({{a, b, 50, .02}}, dt);
	two.SetResistors({{a, b, 100, .01}, {a, b, 100, .01}}, dt);
	for (unsigned int step = 0; step < 1000; ++step)
	{
		one.Step(dt, std::sin(.17 * step));
		two.Step(dt, std::sin(.17 * step));
		Require(std::abs(one.ResistorVoltage(0) - two.ResistorVoltage(0)) < 1e-12,
				"two parallel 100-ohm loads differ from 50 ohms");
	}
	unsigned int rejected = 0;
	try
	{
		one.Step(dt * .9);
	}
	catch (const std::invalid_argument &)
	{
		++rejected;
	}
	try
	{
		two.SetResistors({{a, b, -50, 0}}, dt);
	}
	catch (const std::invalid_argument &)
	{
		++rejected;
	}
	try
	{
		two.VoltagePath(a, {.007, .007, .007});
	}
	catch (const std::invalid_argument &)
	{
		++rejected;
	}
	Require(rejected == 3, "invalid circuit geometry/timestep/value was accepted");
	std::cout << "parallel loads and invalid-input controls: PASS\n";
}
int main(int argc, char **argv)
{
	try
	{
		Integrals();
		Power(argc > 1 && std::string(argv[1]) == "--long" ? 100000 : 5000);
		ParallelAndGuards();
		SoftSourcePower();
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

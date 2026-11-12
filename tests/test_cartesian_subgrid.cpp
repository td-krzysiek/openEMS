/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

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

static double RandomField(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

static void CheckTraceAndStorage()
{
	Grid grid({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu);
	Grid uniform({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu, true);
	Require(grid.CellCount() == 512 + 7 * 64, "refinement does not replace each selected cell with eight children");
	Require(uniform.CellCount() == 4096, "uniform reference dimensions are wrong");
	Require(grid.StorageBytes() < uniform.StorageBytes(), "compact grid has not reduced accounted allocation");
	grid.SetElectric(RandomField);
	const double residual = grid.CurlDivergenceResidual(RandomField);
	std::cout << "D curl residual " << residual << '\n';
	Require(residual < 1e-13, "discrete divergence of curl is not zero");
	double jump = 0;
	for (unsigned int normal = 0; normal < 3; ++normal)
	{
		for (unsigned int side = 0; side < 2; ++side)
		{
			for (unsigned int tangent = 0; tangent < 3; ++tangent)
			{
				if (normal == tangent)
				{
					continue;
				}
				for (unsigned int i = 0; i < 8; ++i)
				{
					Grid::Point a{.00213 + i * .00041, .00327, .00471};
					a[normal] = side ? .006 : .002;
					Grid::Point b = a;
					a[normal] -= 1e-12;
					b[normal] += 1e-12;
					jump = std::max(jump, std::abs(grid.Electric(tangent, a) - grid.Electric(tangent, b)));
				}
			}
		}
	}
	std::cout << "tangential trace max jump " << jump << ", cells " << grid.CellCount() << '/' << uniform.CellCount()
			  << ", accounted bytes " << grid.StorageBytes() << '/' << uniform.StorageBytes() << '\n';
	Require(jump < 1e-7, "tangential electric trace is discontinuous at a box face");
}

static void CheckAnisotropicAndGuards()
{
	Grid grid({6, 8, 10}, {.001, .0014, .0009}, {{1, 2, 3}, {4, 6, 8}}, epsilon, mu);
	grid.SetElectric(RandomField);
	const double dt = .8 * grid.StableTimestep();
	const double initial = grid.Energy(dt);
	for (unsigned int i = 0; i < 2000; ++i)
	{
		grid.Step(dt);
	}
	const double drift = std::abs(grid.Energy(dt) / initial - 1);
	std::cout << "off-centre anisotropic block invariant drift " << drift << '\n';
	Require(drift < 1e-8 && grid.CurlDivergenceResidual(RandomField) < 1e-13,
			"anisotropic/off-centre box conservation failed");
	unsigned int rejected = 0;
	for (unsigned int kind = 0; kind < 4; ++kind)
	{
		try
		{
			if (kind == 0)
			{
				Grid bad({4, 4, 4}, {0, .001, .001}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu);
			}
			else if (kind == 1)
			{
				Grid bad({4, 4, 4}, {.001, .001, .001}, {{0, 1, 1}, {3, 3, 3}}, epsilon, mu);
			}
			else if (kind == 2)
			{
				Grid bad({4, 4, 4}, {.001, .001, .001}, {{1, 1, 1}, {3, 3, 3}}, -epsilon, mu);
			}
			else
			{
				Grid bad({64, 64, 64}, {.001, .001, .001}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu, true);
			}
		}
		catch (const std::invalid_argument &)
		{
			++rejected;
		}
	}
	Require(rejected == 4, "invalid geometry/material/allocation was accepted");
}

static void CheckUniformNative()
{
	Grid compact({4, 4, 4}, {.002, .002, .002}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu, true);
	ContinuousStructure csx;
	CSRectGrid *mesh = csx.GetGrid();
	mesh->SetDeltaUnit(1e-3);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 8; ++i)
		{
			mesh->AddDiscLine(n, i);
		}
	}
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator> op(Operator::New(1));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(op->SetGeometryCSX(&csx), "native mesh failed");
	int boundaries[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(boundaries);
	op->SetTimeStepMethod(1);
	const double dt = .7 * compact.StableTimestep();
	op->SetTimestep(dt);
	Require(op->CalcECOperator() == 0, "native operator failed");
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto field = [](unsigned int n, const Grid::Point &p)
	{ return n == 2 ? std::sin(std::acos(-1.) * p[0] / .008) * std::sin(std::acos(-1.) * p[1] / .008) : 0.; };
	compact.SetElectric(field);
	for (unsigned int x = 0; x < 9; ++x)
	{
		for (unsigned int y = 0; y < 9; ++y)
		{
			for (unsigned int z = 0; z < 8; ++z)
			{
				engine->SetVolt(2, x, y, z, field(2, {x * .001, y * .001, (z + .5) * .001}) * .001);
			}
		}
	}
	for (unsigned int i = 0; i < 200; ++i)
	{
		compact.Step(dt);
		engine->IterateTS(1);
	}
	double peakE = 0, errorE = 0, peakH = 0, errorH = 0;
	for (unsigned int x = 1; x < 8; ++x)
	{
		for (unsigned int y = 1; y < 8; ++y)
		{
			for (unsigned int z = 1; z < 8; ++z)
			{
				const unsigned int pos[3] = {x, y, z};
				for (unsigned int n = 0; n < 3; ++n)
				{
					double a[3], b[3];
					unsigned int p[3] = {x, y, z};
					op->GetYeeCoords(n, p, a, false);
					op->GetYeeCoords(n, p, b, true);
					const Grid::Point e{a[0] * .001, a[1] * .001, a[2] * .001};
					const Grid::Point h{b[0] * .001, b[1] * .001, b[2] * .001};
					const double nativeE = engine->GetVolt(n, x, y, z) / op->GetEdgeLength(n, pos);
					const double nativeH = engine->GetCurr(n, x, y, z) / op->GetEdgeLength(n, pos, true);
					peakE = std::max(peakE, std::abs(nativeE));
					peakH = std::max(peakH, std::abs(nativeH));
					errorE = std::max(errorE, std::abs(nativeE - compact.Electric(n, e)));
					errorH = std::max(errorH, std::abs(nativeH - compact.Magnetic(n, h)));
				}
			}
		}
	}
	std::cout << "uniform vs native: E " << errorE / peakE << ", H " << errorH / peakH << '\n';
	Require(peakE > 0 && peakH > 0 && errorE / peakE < 2e-4 && errorH / peakH < 2e-4,
			"uniform core does not recover native Yee propagation");
}

static void CheckEnergy(unsigned int steps)
{
	Grid grid({4, 4, 4}, {.002, .002, .002}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu);
	grid.SetElectric(RandomField);
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return RandomField(n, p) / 376.730313668; });
	const double dt = .8 * grid.StableTimestep();
	const double initial = grid.Energy(dt);
	double error = 0;
	for (unsigned int step = 0; step < steps; ++step)
	{
		grid.Step(dt);
		if (step % 100 == 0)
		{
			error = std::max(error, std::abs(grid.Energy(dt) / initial - 1));
			Require(std::isfinite(error) && error < 1e-8, "discrete energy invariant is not conserved");
		}
	}
	std::cout << steps << " steps: relative invariant drift " << error << '\n';
	bool rejected = false;
	try
	{
		grid.Step(grid.StableTimestep() * 1.1);
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "unsafe step accepted");
}

static double CavityFrequency(unsigned int cells, unsigned int component)
{
	constexpr double length = .012;
	const double pi = std::acos(-1.);
	Grid grid({cells, cells, cells}, {length / cells, length / cells, length / cells},
			  {{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, epsilon, mu);
	grid.SetElectric(
		[&](unsigned int n, const Grid::Point &p)
		{
			return n == component ? std::sin(pi * p[(n + 1) % 3] / length) * std::sin(pi * p[(n + 2) % 3] / length)
								  : 0.;
		});
	const double dt = .8 * grid.StableTimestep();
	const Grid::Point probe{.00317, .00431, .00547};
	double previous = grid.Electric(component, probe);
	std::vector<double> crossing;
	for (unsigned int step = 1; step <= 4000 * cells / 6; ++step)
	{
		grid.Step(dt);
		const double value = grid.Electric(component, probe);
		if (previous <= 0 && value > 0)
		{
			crossing.push_back(step - value / (value - previous));
		}
		previous = value;
	}
	Require(crossing.size() > 15, "no cavity mode recorded");
	const double frequency = (crossing.size() - 1) / ((crossing.back() - crossing.front()) * dt);
	const double exact = std::sqrt(2.) / (2 * length * std::sqrt(epsilon * mu));
	std::cout << "cavity base " << cells << ", E" << component << ": frequency " << frequency
			  << ", relative analytic error " << frequency / exact - 1 << '\n';
	return frequency / exact - 1;
}

static double PulseError(unsigned int cells, unsigned int propagation)
{
	constexpr double length = .012;
	const Grid::Index dimensions{cells, cells, cells};
	const Grid::Point spacing{length / cells, length / cells, length / cells};
	const Grid::Box box{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}};
	Grid grid(dimensions, spacing, box, epsilon, mu);
	Grid reference(dimensions, spacing, box, epsilon, mu, true);
	const unsigned int polarization = (propagation + 1) % 3;
	auto seed = [&](unsigned int n, const Grid::Point &p)
	{
		return n == polarization ? std::exp(-std::pow((p[propagation] - .003) / .002, 2)) *
									   std::sin(std::acos(-1.) * p[propagation] / length) *
									   std::sin(std::acos(-1.) * p[(propagation + 2) % 3] / length)
								 : 0.;
	};
	grid.SetElectric(seed);
	reference.SetElectric(seed);
	const double dt = .8 * std::min(grid.StableTimestep(), reference.StableTimestep());
	Grid::Point point{.006, .006, .006};
	point[propagation] = .009;
	double error = 0, peak = 0;
	const unsigned int steps = unsigned(1.1 * length * std::sqrt(epsilon * mu) / dt);
	for (unsigned int step = 0; step < steps; ++step)
	{
		grid.Step(dt);
		reference.Step(dt);
		const double a = reference.Electric(polarization, point);
		peak = std::max(peak, std::abs(a));
		error = std::max(error, std::abs(a - grid.Electric(polarization, point)));
	}
	std::cout << "pulse axis " << propagation << ", base " << cells << ": peak relative error " << error / peak << '\n';
	Require(peak > 1e-3, "no pulse at reference receiver");
	return error / peak;
}

int main(int argc, char **argv)
{
	try
	{
		if (argc == 2 && std::string(argv[1]) == "--mesh")
		{
			Grid grid({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu);
			std::cout << "{\"baseCells\":[8,8,8],\"spacingMm\":[1,1,1],\"refinedBoxMm\":[2,2,2,6,6,6],\"leaves\":[";
			bool first = true;
			for (const auto &cell : grid.Leaves())
			{
				if (!first)
				{
					std::cout << ',';
				}
				first = false;
				std::cout << '[' << cell.lo[0] / 2. << ',' << cell.lo[1] / 2. << ',' << cell.lo[2] / 2. << ','
						  << cell.size / 2. << ']';
			}
			std::cout << "],\"leafCount\":" << grid.CellCount() << ",\"accountedBytes\":" << grid.StorageBytes()
					  << "}\n";
			return 0;
		}
		CheckTraceAndStorage();
		CheckAnisotropicAndGuards();
		CheckUniformNative();
		CheckEnergy(argc == 2 && std::string(argv[1]) == "--long" ? 1000000 : 20000);
		const double coarse = CavityFrequency(6, 2);
		const double fine = CavityFrequency(12, 2);
		Require(std::abs(fine) < .005 && std::abs(fine) < .75 * std::abs(coarse), "3-D cavity is not converging");
		Require(std::abs(CavityFrequency(6, 0) - coarse) < 2e-3 && std::abs(CavityFrequency(6, 1) - coarse) < 2e-3,
				"coordinate permutations produce inconsistent cavity frequencies");
		for (unsigned int direction = 0; direction < 3; ++direction)
		{
			const double a = PulseError(6, direction);
			const double b = PulseError(12, direction);
			const double c = PulseError(24, direction);
			Require(b < a && c < .75 * b && c < .1,
					"transmitted pulse is not converging to the uniform fine reference");
		}
		Require(PulseError(36, 0) < .05, "well-resolved pulse differs by more than five percent from uniform fine");
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	return 0;
}

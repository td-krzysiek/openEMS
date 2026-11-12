/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "CSPrimBox.h"
#include "CSPropMaterial.h"
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

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

static Grid::Axes Graded(Grid::Index cells, Grid::Point origin = {0, 0, 0}, Grid::Point lengths = {.012, .008, .008})
{
	const Grid::Point warp{.18, -.15, .12};
	Grid::Axes lines;
	const double pi = std::acos(-1.);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= cells[n]; ++i)
		{
			const double t = double(i) / cells[n];
			lines[n].push_back(origin[n] + lengths[n] * (t + warp[n] * std::sin(2 * pi * t) / (2 * pi)));
		}
	}
	return lines;
}

static std::vector<Grid::Region> Regions()
{
	return {{{{1, 2, 2}, {4, 6, 6}}, {2, 2, 4}}, {{{7, 2, 2}, {10, 6, 6}}, {4, 1, 2}}};
}

static double Seed(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

static void CoordinatesAndTranslation(unsigned int steps)
{
	const Grid::Point origin{-.004, -.002, .001}, shift{.031, -.023, .017};
	const auto lines = Graded({12, 8, 8}, origin);
	Grid grid(lines, epsilon, mu, Regions());
	Require(grid.BaseLines() == lines && grid.CellCount() == 1824,
			"nonuniform coordinates or local cell count changed");
	double volume = 0;
	const auto scale = grid.LatticeRatio();
	for (std::size_t i = 0; i < grid.CellCount(); ++i)
	{
		const auto &leaf = grid.Leaves()[i];
		const auto start = grid.LeafStart(i), length = grid.LeafLengths(i);
		double local = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			const unsigned int root = unsigned(leaf.lo[n]) / scale[n];
			const double base = lines[n][root + 1] - lines[n][root];
			const double wantStart = lines[n][root] + base * double(leaf.lo[n] % int(scale[n])) / scale[n];
			Require(std::abs(start[n] - wantStart) < 1e-16 &&
						std::abs(length[n] - base * leaf.size[n] / scale[n]) < 1e-16,
					"local physical dimensions do not follow base interval widths");
			Require(start[n] >= lines[n][root] - 1e-16 && start[n] + length[n] <= lines[n][root + 1] + 1e-16,
					"refined leaf escapes its own base interval");
			local *= length[n];
		}
		volume += local;
	}
	Require(std::abs(volume / (.012 * .008 * .008) - 1) < 1e-12, "nonuniform leaves do not tile the domain");
	Grid::Axes translated = lines;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (double &v : translated[n])
		{
			v += shift[n];
		}
	}
	Grid other(translated, epsilon, mu, Regions());
	auto seed = [&](unsigned int n, Grid::Point p)
	{
		for (unsigned int k = 0; k < 3; ++k)
		{
			p[k] -= origin[k];
		}
		return Seed(n, p);
	};
	grid.SetElectric(seed);
	grid.SetMagnetic([&](unsigned int n, const Grid::Point &p) { return seed(n, p) / 376.730313668; });
	other.SetElectric(
		[&](unsigned int n, Grid::Point p)
		{
			for (unsigned int k = 0; k < 3; ++k)
			{
				p[k] -= shift[k];
			}
			return seed(n, p);
		});
	other.SetMagnetic(
		[&](unsigned int n, Grid::Point p)
		{
			for (unsigned int k = 0; k < 3; ++k)
			{
				p[k] -= shift[k];
			}
			return seed(n, p) / 376.730313668;
		});
	double jump = 0;
	for (const auto &region : grid.Regions())
	{
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
					Grid::Point a;
					for (unsigned int n = 0; n < 3; ++n)
					{
						a[n] = .63 * lines[n][region.box.start[n]] + .37 * lines[n][region.box.stop[n]];
					}
					a[normal] = lines[normal][side ? region.box.stop[normal] : region.box.start[normal]];
					auto b = a;
					a[normal] -= 1e-12;
					b[normal] += 1e-12;
					jump = std::max(jump, std::abs(grid.Electric(tangent, a) - grid.Electric(tangent, b)));
				}
			}
		}
	}
	const double divcurl = grid.CurlDivergenceResidual(seed);
	Require(jump < 1e-7 && divcurl < 1e-13, "nonuniform interface continuity or div(curl) failed");
	const double dt = .8 * std::min(grid.StableTimestep(), other.StableTimestep());
	const double energy = grid.Energy(dt);
	double drift = 0, errorE = 0, errorH = 0, peakE = 0, peakH = 0;
	for (unsigned int i = 0; i < steps; ++i)
	{
		grid.Step(dt);
		other.Step(dt);
		if (i % 100 == 0)
		{
			drift = std::max(drift, std::abs(grid.Energy(dt) / energy - 1));
			Require(std::isfinite(drift) && drift < 1e-8, "nonuniform invariant drift exceeds limit");
			for (const Grid::Point relative :
				 {Grid::Point{.00213, .00317, .00471}, {.00513, .00271, .00367}, {.00811, .00513, .00227}})
			{
				Grid::Point a, b;
				for (unsigned int n = 0; n < 3; ++n)
				{
					a[n] = origin[n] + relative[n];
					b[n] = a[n] + shift[n];
				}
				for (unsigned int n = 0; n < 3; ++n)
				{
					const double e = grid.Electric(n, a), h = grid.Magnetic(n, a);
					peakE = std::max(peakE, std::abs(e));
					peakH = std::max(peakH, std::abs(h));
					errorE = std::max(errorE, std::abs(e - other.Electric(n, b)));
					errorH = std::max(errorH, std::abs(h - other.Magnetic(n, b)));
				}
			}
		}
	}
	Require(peakE > 0 && peakH > 0 && errorE / peakE < 1e-8 && errorH / peakH < 1e-8,
			"coordinate translation changed fields");
	std::cout << "graded multi-region: " << steps << " steps, invariant " << drift << ", trace " << jump << ", divcurl "
			  << divcurl << ", translation E/H " << errorE / peakE << '/' << errorH / peakH << '\n';
}

static void NativeGradedYee(bool dielectric = false)
{
	const auto lines = Graded({6, 6, 6}, {.001, -.002, .003}, {.008, .008, .008});
	const Grid::Index ratio{2, 1, 4};
	const std::vector<Grid::Region> regions{{{{1, 1, 1}, {5, 5, 5}}, ratio}};
	const double boundary = lines[0][3];
	auto medium = [&](const Grid::Point &start, const Grid::Point &stop)
	{
		return Grid::Medium{epsilon * (dielectric && (start[0] + stop[0]) / 2 > boundary ? 3.556 : 1.), mu};
	};
	Grid compact(lines, epsilon, mu, regions, true, medium);
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(1.);
	if (dielectric)
	{
		auto *material = new CSPropMaterial(csx.GetParameterSet());
		material->SetEpsilon(3.556);
		csx.AddProperty(material);
		auto *box = new CSPrimBox(csx.GetParameterSet(), material);
		for (unsigned int n = 0; n < 3; ++n)
		{
			box->SetCoord(2 * n, n == 0 ? boundary : lines[n].front());
			box->SetCoord(2 * n + 1, lines[n].back());
		}
	}
	Grid::Axes fine;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i < 6; ++i)
		{
			for (unsigned int k = 0; k < ratio[n]; ++k)
			{
				fine[n].push_back(lines[n][i] + (lines[n][i + 1] - lines[n][i]) * k / ratio[n]);
			}
		}
		fine[n].push_back(lines[n].back());
		for (double v : fine[n])
		{
			csx.GetGrid()->AddDiscLine(n, v);
		}
	}
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator> op(Operator::New(1));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(op->SetGeometryCSX(&csx), "native graded mesh rejected");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	const double dt = .7 * compact.StableTimestep();
	op->SetTimestep(dt);
	Require(op->CalcECOperator() == 0, "native graded operator failed");
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto seed = [&](unsigned int n, const Grid::Point &p)
	{
		const double pi = std::acos(-1.);
		return std::sin(pi * (p[(n + 1) % 3] - lines[(n + 1) % 3].front()) / .008) *
			   std::sin(pi * (p[(n + 2) % 3] - lines[(n + 2) % 3].front()) / .008) *
			   std::cos(2 * pi * (p[n] - lines[n].front()) / .008);
	};
	compact.SetElectric(seed);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int x = 0; x < fine[0].size(); ++x)
		{
			for (unsigned int y = 0; y < fine[1].size(); ++y)
			{
				for (unsigned int z = 0; z < fine[2].size(); ++z)
				{
					unsigned int pos[3] = {x, y, z};
					if (pos[n] + 1 == fine[n].size())
					{
						continue;
					}
					double p[3];
					op->GetYeeCoords(n, pos, p, false);
					engine->SetVolt(n, pos, seed(n, {p[0], p[1], p[2]}) * op->GetEdgeLength(n, pos));
				}
			}
		}
	}
	for (unsigned int i = 0; i < 200; ++i)
	{
		compact.Step(dt);
		engine->IterateTS(1);
	}
	double errorE = 0, errorH = 0, peakE = 0, peakH = 0;
	for (unsigned int x = 1; x + 1 < fine[0].size(); ++x)
	{
		for (unsigned int y = 1; y + 1 < fine[1].size(); ++y)
		{
			for (unsigned int z = 1; z + 1 < fine[2].size(); ++z)
			{
				unsigned int pos[3] = {x, y, z};
				for (unsigned int n = 0; n < 3; ++n)
				{
					double e[3], h[3];
					op->GetYeeCoords(n, pos, e, false);
					op->GetYeeCoords(n, pos, h, true);
					const double nativeE = engine->GetVolt(n, pos) / op->GetEdgeLength(n, pos);
					const double nativeH = engine->GetCurr(n, pos) / op->GetEdgeLength(n, pos, true);
					peakE = std::max(peakE, std::abs(nativeE));
					peakH = std::max(peakH, std::abs(nativeH));
					errorE = std::max(errorE, std::abs(nativeE - compact.Electric(n, {e[0], e[1], e[2]})));
					errorH = std::max(errorH, std::abs(nativeH - compact.Magnetic(n, {h[0], h[1], h[2]})));
				}
			}
		}
	}
	std::cout << "graded uniform-fine " << (dielectric ? "dielectric " : "vacuum ") << "vs native Yee E/H "
			  << errorE / peakE << '/' << errorH / peakH << '\n';
	Require(peakE > 0 && peakH > 0 && errorE / peakE < 2e-4 && errorH / peakH < 2e-4,
			"graded uniform-fine core differs from native Yee");
}

static void InvalidAxes()
{
	const auto valid = Graded({6, 6, 6});
	for (unsigned int kind = 0; kind < 6; ++kind)
	{
		auto lines = valid;
		if (kind == 0)
		{
			lines[0][2] = lines[0][1];
		}
		if (kind == 1)
		{
			lines[1][2] = lines[1][1] - .001;
		}
		if (kind == 2)
		{
			lines[2][3] = std::numeric_limits<double>::infinity();
		}
		if (kind == 3)
		{
			lines[2][3] = std::numeric_limits<double>::quiet_NaN();
		}
		if (kind == 4)
		{
			lines[0].resize(4);
		}
		if (kind == 5)
		{
			for (unsigned int i = 0; i < 7; ++i)
			{
				lines[0][i] = 1e16 + 2. * i;
			}
		}
		bool rejected = false;
		try
		{
			Grid bad(lines, epsilon, mu, {{{{1, 1, 1}, {5, 5, 5}}, {2, 2, 2}}});
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected, "invalid or unrepresentable nonuniform mesh accepted");
	}
	// A thin but representable interval is valid and must reduce the step.
	auto thin = valid;
	thin[2][2] = thin[2][1] + 1e-6;
	Grid normal(valid, epsilon, mu, {{{{1, 1, 1}, {5, 5, 5}}, {2, 2, 2}}});
	Grid narrow(thin, epsilon, mu, {{{{1, 1, 1}, {5, 5, 5}}, {2, 2, 2}}});
	Require(narrow.StableTimestep() < normal.StableTimestep() / 10, "thin interval did not constrain the coupled step");
	narrow.SetElectric(Seed);
	const double dt = .8 * narrow.StableTimestep(), initial = narrow.Energy(dt);
	for (unsigned int i = 0; i < 1000; ++i)
	{
		narrow.Step(dt);
	}
	Require(std::abs(narrow.Energy(dt) / initial - 1) < 1e-8, "thin interval is unstable at the coupled timestep");
	std::cout << "thin interval timestep ratio " << narrow.StableTimestep() / normal.StableTimestep() << '\n';
}

static double Packet(Grid::Index cells)
{
	const auto axes = Graded(cells, {-.002, .001, -.003}, {.012, .012, .012});
	const std::vector<Grid::Region> regions{
		{{{cells[0] / 6, cells[1] / 3, cells[2] / 3}, {cells[0] / 3, 2 * cells[1] / 3, 2 * cells[2] / 3}}, {2, 2, 2}},
		{{{2 * cells[0] / 3, cells[1] / 3, cells[2] / 3}, {5 * cells[0] / 6, 2 * cells[1] / 3, 2 * cells[2] / 3}},
		 {4, 2, 1}}};
	Grid grid(axes, epsilon, mu, regions), reference(axes, epsilon, mu, regions, true);
	auto seed = [&](unsigned int n, const Grid::Point &p)
	{
		const double x = p[0] - axes[0].front(), z = p[2] - axes[2].front();
		return n == 1 ? std::exp(-std::pow((x - .003) / .002, 2)) * std::sin(std::acos(-1.) * x / .012) *
							std::sin(std::acos(-1.) * z / .012)
					  : 0.;
	};
	grid.SetElectric(seed);
	reference.SetElectric(seed);
	const double dt = .8 * std::min(grid.StableTimestep(), reference.StableTimestep());
	const unsigned int steps = unsigned(1.1 * .012 * std::sqrt(epsilon * mu) / dt);
	const Grid::Point receiver{axes[0].front() + .011, axes[1].front() + .006, axes[2].front() + .006};
	double error = 0, peak = 0;
	for (unsigned int i = 0; i < steps; ++i)
	{
		grid.Step(dt);
		reference.Step(dt);
		const double value = reference.Electric(1, receiver);
		peak = std::max(peak, std::abs(value));
		error = std::max(error, std::abs(value - grid.Electric(1, receiver)));
	}
	Require(peak > 1e-3, "no packet beyond graded refinement regions");
	std::cout << "graded packet base " << cells[0] << 'x' << cells[1] << 'x' << cells[2] << " error " << error / peak
			  << '\n';
	return error / peak;
}

int main(int argc, char **argv)
{
	try
	{
		CoordinatesAndTranslation(argc == 2 && std::string(argv[1]) == "--long" ? 100000 : 4000);
		NativeGradedYee();
		NativeGradedYee(true);
		InvalidAxes();
		const double a = Packet({6, 6, 6}), b = Packet({12, 12, 12}), c = Packet({24, 24, 24}),
					 d = Packet({36, 36, 36});
		const double resolved = Packet({42, 36, 36});
		Require(b < a && c < b && d < c && resolved < d && resolved < .05,
				"graded multi-region packet does not converge within five percent");
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	return 0;
}

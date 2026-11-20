/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "ContinuousStructure.h"
#include "CSPrimBox.h"
#include "CSPrimCylinder.h"
#include "CSPrimCylindricalShell.h"
#include "CSPropMetal.h"

#include <algorithm>
#include <cmath>
#include <iostream>
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

static Grid::Axes Lines(unsigned int count)
{
	Grid::Axes axes;
	for (auto &axis : axes)
	{
		for (unsigned int i = 0; i <= count; ++i)
		{
			axis.push_back(.012 * i / count);
		}
	}
	return axes;
}

static double Seed(unsigned int n, const Grid::Point &p)
{
	return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n) +
		   .3 * std::cos(4381 * p[0] - 2271 * p[1] + 411 * p[2] + 2 * n);
}

static Grid::PECSelector BoxPEC(Grid::Point lo, Grid::Point hi)
{
	return [=](unsigned int axis, const Grid::Point &start, const Grid::Point &stop)
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double p = (start[n] + stop[n]) / 2;
			if (p < lo[n] - 1e-15 || p > hi[n] + 1e-15 || (lo[n] == hi[n] && n == axis))
			{
				return false;
			}
		}
		return true;
	};
}

static void ClosedCavityAndInvariant(unsigned int steps)
{
	const auto axes = Lines(6);
	const std::vector<Grid::Region> regions{{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}};
	const auto sheet = BoxPEC({.006, 0, 0}, {.006, .012, .012});
	Grid grid(axes, epsilon, mu, regions, false, {}, sheet);
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
	Require(std::isfinite(drift) && drift < 1e-8, "PEC elimination broke lossless invariant");
	Require(grid.CurlDivergenceResidual(Seed) < 1e-13, "PEC elimination broke div(curl)");
	Require(grid.PECTraceCount() > 0, "no PEC traces classified");
	Require(std::abs(grid.Electric(1, {.006, .00517, .00471})) < 1e-12 &&
				std::abs(grid.Electric(2, {.006, .00517, .00471})) < 1e-12 &&
				std::abs(grid.Magnetic(0, {.006, .00517, .00471})) < 1e-12,
			"PEC tangential E or normal B is nonzero");
	std::cout << "PEC sheet crossing refinement: invariant " << steps << " steps " << drift << '\n';
}

static void SolidInvariant(unsigned int steps)
{
	const auto axes = Lines(6);
	const std::vector<Grid::Region> regions{{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}};
	for (bool crossing : {false, true})
	{
		const Grid::Point lo = crossing ? Grid::Point{.002, .004, .004} : Grid::Point{.005, .005, .005};
		const Grid::Point hi = crossing ? Grid::Point{.006, .008, .008} : Grid::Point{.007, .007, .007};
		Grid grid(axes, epsilon, mu, regions, false, {}, BoxPEC(lo, hi));
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
		const Grid::Point inside =
			crossing ? Grid::Point{.00517, .00531, .00547} : Grid::Point{.00617, .00631, .00647};
		for (unsigned int n = 0; n < 3; ++n)
		{
			Require(std::abs(grid.Electric(n, inside)) < 1e-12 && std::abs(grid.Magnetic(n, inside)) < 1e-12,
					"solid PEC stores internal fields");
		}
		Require(std::isfinite(drift) && drift < 1e-8, "solid PEC invariant failed");
		Require(grid.CurlDivergenceResidual(Seed) < 1e-13, "solid PEC broke div(curl)");
		std::cout << "solid PEC " << (crossing ? "crossing" : "inside") << " invariant " << steps << " steps "
				  << drift << '\n';
	}
}

static double Cavity(unsigned int cells)
{
	const double pi = std::acos(-1.), length = .012;
	Grid grid(
		Lines(cells), epsilon, mu,
		{{{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, {2, 2, 2}}},
		false, {}, BoxPEC({length / 2, 0, 0}, {length / 2, length, length}));
	grid.SetElectric(
		[&](unsigned int n, const Grid::Point &p)
		{
			return n == 2 && p[0] < length / 2
					   ? std::sin(2 * pi * p[0] / length) * std::sin(pi * p[1] / length)
					   : 0.;
		});
	const double dt = .8 * grid.StableTimestep();
	const Grid::Point receiver{.00317, .00431, .00547};
	double previous = grid.Electric(2, receiver), leakage = 0;
	std::vector<double> crossings;
	for (unsigned int step = 1; step <= 4000 * cells / 6; ++step)
	{
		grid.Step(dt);
		const double value = grid.Electric(2, receiver);
		if (previous <= 0 && value > 0)
		{
			crossings.push_back(step - value / (value - previous));
		}
		previous = value;
		leakage = std::max(leakage, std::abs(grid.Electric(2, {.00917, .00431, .00547})));
	}
	Require(crossings.size() > 15, "no split-cavity mode observed");
	Require(leakage < 1e-12, "PEC wall transmitted fields into disconnected cavity");
	const double frequency = (crossings.size() - 1) / ((crossings.back() - crossings.front()) * dt);
	const double exact = std::sqrt(5.) / (2 * length * std::sqrt(epsilon * mu));
	const double error = frequency / exact - 1;
	std::cout << "split PEC cavity " << cells << " frequency error " << error << " leakage " << leakage
			  << '\n';
	return std::abs(error);
}

static void UnsupportedInterface()
{
	const auto axes = Lines(6);
	const std::vector<Grid::Region> regions{{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}};
	for (bool opening : {false, true})
	{
		bool rejected = false;
		try
		{
			auto mask = [=](unsigned int n, const Grid::Point &a, const Grid::Point &b)
			{
				const double y = (a[1] + b[1]) / 2, z = (a[2] + b[2]) / 2;
				const bool plane = std::abs(a[0] - .004) < 1e-15 && std::abs(b[0] - .004) < 1e-15;
				if (!plane || n == 0)
				{
					return false;
				}
				// A fine strip, or a fine aperture inside an otherwise PEC face,
				// which cannot be represented by the adjacent coarse trace.
				const bool strip = y > .0044 && y < .0056 && z > .0044 && z < .0056;
				return opening ? !strip : strip;
			};
			Grid bad(axes, epsilon, mu, regions, false, {}, mask);
		}
		catch (const std::invalid_argument &error)
		{
			rejected = std::string(error.what()).find("coarse interface trace") != std::string::npos;
		}
		Require(rejected, "unrepresentable PEC strip/opening accepted at the interface");
	}
}

static CSPrimBox *AddBox(ContinuousStructure &csx, CSProperties *property, Grid::Point lo, Grid::Point hi)
{
	auto *box = new CSPrimBox(csx.GetParameterSet(), property);
	box->SetPriority(10);
	for (unsigned int n = 0; n < 3; ++n)
	{
		box->SetCoord(2 * n, lo[n]);
		box->SetCoord(2 * n + 1, hi[n]);
	}
	return box;
}

static void NativeYee(unsigned int shape)
{
	const auto axes = Lines(6);
	const double length = .012;
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(1.);
	Grid::Axes fine;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i < 6; ++i)
		{
			fine[n].push_back(axes[n][i]);
			fine[n].push_back((axes[n][i] + axes[n][i + 1]) / 2);
		}
		fine[n].push_back(length);
		for (double v : fine[n])
		{
			csx.GetGrid()->AddDiscLine(n, v);
		}
	}
	auto *metal = new CSPropMetal(csx.GetParameterSet());
	csx.AddProperty(metal);
	Grid::PECSelector mask;
	if (shape < 2)
	{
		const Grid::Point lo = shape == 0 ? Grid::Point{.004, .004, .004} : Grid::Point{.006, 0, 0};
		const Grid::Point hi = shape == 0 ? Grid::Point{.008, .008, .006} : Grid::Point{.006, length, length};
		AddBox(csx, metal, lo, hi);
		mask = BoxPEC(lo, hi);
	}
	else
	{
		CSPrimCylinder *cylinder;
		if (shape == 2)
		{
			cylinder = new CSPrimCylinder(csx.GetParameterSet(), metal);
		}
		else
		{
			auto *shell = new CSPrimCylindricalShell(csx.GetParameterSet(), metal);
			shell->SetShellWidth(.001);
			cylinder = shell;
		}
		cylinder->SetRadius(.0025);
		for (unsigned int n = 0; n < 3; ++n)
		{
			cylinder->SetCoord(2 * n, n == 2 ? .002 : .006);
			cylinder->SetCoord(2 * n + 1, n == 2 ? .010 : .006);
		}
		mask = [=](unsigned int, const Grid::Point &a, const Grid::Point &b)
		{
			const Grid::Point p{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
			const double radius = std::hypot(p[0] - .006, p[1] - .006);
			return p[2] >= .002 - 1e-15 && p[2] <= .010 + 1e-15 &&
				   (shape == 2 ? radius <= .0025 + 1e-15 : radius >= .002 - 1e-15 && radius <= .003 + 1e-15);
		};
	}
	Grid compact(axes, epsilon, mu, {{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}}, true, {}, mask);
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator> op(Operator::New(1));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(csx.Update().empty(), "native PEC fixture geometry update failed");
	Require(op->SetGeometryCSX(&csx), "native PEC geometry rejected");
	int bc[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(bc);
	const double dt = .7 * compact.StableTimestep();
	op->SetTimestep(dt);
	Require(op->CalcECOperator() == 0, "native PEC operator failed");
	std::unique_ptr<Engine> engine(op->CreateEngine());
	auto seed = [&](unsigned int n, const Grid::Point &p)
	{
		return std::sin(std::acos(-1.) * p[(n + 1) % 3] / length) *
			   std::sin(std::acos(-1.) * p[(n + 2) % 3] / length) *
			   std::cos(2 * std::acos(-1.) * p[n] / length);
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
					Grid::Point a{p[0], p[1], p[2]}, b = a;
					a[n] = fine[n][pos[n]];
					b[n] = fine[n][pos[n] + 1];
					engine->SetVolt(
						n, pos, mask(n, a, b) ? 0 : seed(n, {p[0], p[1], p[2]}) * op->GetEdgeLength(n, pos));
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
	std::cout << "PEC shape " << shape << " native Yee E/H " << errorE / peakE << '/' << errorH / peakH
			  << '\n';
	Require(peakE > 0 && peakH > 0 && errorE / peakE < 2e-4 && errorH / peakH < 2e-4,
			"PEC core differs from native Yee");
}

static void NativeContact(bool gap)
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(1.);
	Grid::Axes axes;
	for (unsigned int n = 0; n < 2; ++n)
	{
		for (unsigned int i = 0; i <= 8; ++i)
		{
			axes[n].push_back(.001 * i);
		}
	}
	axes[2] = {0, .001, .002, .002018, .003, .003982, .004, .004018, .005, .006, .007, .008};
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (double value : axes[n])
		{
			csx.GetGrid()->AddDiscLine(n, value);
		}
	}
	auto *metal = new CSPropMetal(csx.GetParameterSet());
	csx.AddProperty(metal);
	AddBox(csx, metal, {.001, .001, .002}, {.007, .007, .002018});
	AddBox(csx, metal, {.003, .003, .004}, {.005, .005, .004018});
	auto *via = new CSPrimCylindricalShell(csx.GetParameterSet(), metal);
	via->SetPriority(20);
	via->SetRadius(.00025);
	via->SetShellWidth(.0002);
	const double end = gap ? .003982 : .004;
	for (unsigned int n = 0; n < 3; ++n)
	{
		via->SetCoord(2 * n, n == 2 ? .002018 : .004);
		via->SetCoord(2 * n + 1, n == 2 ? end : .004);
	}
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(
		Operator_CartesianSubgrid::New({.002, .002, .001, .006, .006, .006}, {4, 4, 1}));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(csx.Update().empty(), "native Cu18 fixture geometry update failed");
	Require(op->SetGeometryCSX(&csx), "native Cu18 contact mesh rejected");
	Require(op->CalcECOperator() == 0, "native Cu18 contact operator failed");
	auto &grid = op->Grid();
	std::size_t checked = 0, blocked = 0;
	for (std::size_t cell = 0; cell < grid.CellCount(); ++cell)
	{
		const auto lo = grid.LeafStart(cell), size = grid.LeafLengths(cell);
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int a = 0; a < 2; ++a)
			{
				for (unsigned int b = 0; b < 2; ++b)
				{
					Grid::Point point = lo;
					point[n] += .5 * size[n];
					point[(n + 1) % 3] += a * size[(n + 1) % 3];
					point[(n + 2) % 3] += b * size[(n + 2) % 3];
					auto box = [&](Grid::Point first, Grid::Point last)
					{
						for (unsigned int k = 0; k < 3; ++k)
						{
							if (point[k] < first[k] - 1e-15 || point[k] > last[k] + 1e-15)
							{
								return false;
							}
						}
						return true;
					};
					const double radius = std::hypot(point[0] - .004, point[1] - .004);
					const bool wall = point[2] >= .002018 - 1e-15 && point[2] <= end + 1e-15 &&
									  radius >= .00015 - 1e-15 && radius <= .00035 + 1e-15;
					const bool expected = box({.001, .001, .002}, {.007, .007, .002018}) ||
										  box({.003, .003, .004}, {.005, .005, .004018}) || wall;
					if (grid.LeafEdgeIsPEC(cell, n, a, b) != expected)
					{
						std::cerr << "occupancy mismatch at " << point[0] << ',' << point[1] << ','
								  << point[2] << " axis " << n << " expected " << expected << " actual "
								  << grid.LeafEdgeIsPEC(cell, n, a, b) << '\n';
						throw std::runtime_error("native Cu18/barrel occupancy mismatch");
					}
					++checked;
					blocked += expected;
				}
			}
		}
	}
	grid.SetElectric(
		[](unsigned int n, const Grid::Point &)
		{
			return n == 2 ? 1. : 0.;
		});
	const Grid::Point contact{.00425, .004, .003991}, bore{.004, .004, .003991};
	const double value = grid.Electric(2, contact);
	Require(std::abs(value - (gap ? 1 : 0)) < 1e-12,
			"18um contact gap was shorted or a connected barrel opened");
	Require(std::abs(grid.Electric(2, bore) - 1) < 1e-12, "hollow barrel was replaced by a solid cylinder");
	std::cout << "native Cu18/18 " << (gap ? "18um GAP" : "CONNECTED") << ": " << checked << " edge checks, "
			  << blocked << " blocked samples; wall/bore E " << value << '/' << grid.Electric(2, bore)
			  << '\n';
}

int main(int argc, char **argv)
{
	try
	{
		NativeContact(false);
		NativeContact(true);
		SolidInvariant(argc == 2 && std::string(argv[1]) == "--long" ? 100000 : 4000);
		UnsupportedInterface();
		ClosedCavityAndInvariant(argc == 2 && std::string(argv[1]) == "--long" ? 100000 : 4000);
		const double coarse = Cavity(6), fine = Cavity(12);
		Require(fine < .01 && fine < .8 * coarse, "split PEC cavity does not converge");
		for (unsigned int shape = 0; shape < 4; ++shape)
		{
			NativeYee(shape);
		}
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
	return 0;
}

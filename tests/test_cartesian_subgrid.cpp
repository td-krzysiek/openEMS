/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "ContinuousStructure.h"
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "FDTD/operator_cartesiansubgrid.h"

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

static void CheckTraceAndStorage(Grid::Index ratio = {2, 2, 2})
{
	Grid grid({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu, false, ratio);
	Grid uniform({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu, true, ratio);
	const unsigned int children = ratio[0] * ratio[1] * ratio[2];
	Require(grid.CellCount() == 512 + (children - 1) * 64, "refinement leaf count does not match axis ratios");
	Require(uniform.CellCount() == 512 * children, "uniform reference dimensions are wrong");
	Require(grid.StorageBytes() <= uniform.StorageBytes(), "compact grid has not reduced accounted allocation");
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
	std::cout << "ratio " << ratio[0] << 'x' << ratio[1] << 'x' << ratio[2] << ": tangential trace max jump " << jump
			  << ", cells " << grid.CellCount() << '/' << uniform.CellCount() << ", accounted bytes "
			  << grid.StorageBytes() << '/' << uniform.StorageBytes() << '\n';
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

static void CheckUniformNative(Grid::Index ratio = {2, 2, 2})
{
	Grid compact({4, 4, 4}, {.002, .002, .002}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu, true, ratio);
	ContinuousStructure csx;
	CSRectGrid *mesh = csx.GetGrid();
	mesh->SetDeltaUnit(1e-3);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 4 * ratio[n]; ++i)
		{
			mesh->AddDiscLine(n, 2. * i / ratio[n]);
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
	{
		const double pi = std::acos(-1.);
		return std::sin(pi * p[(n + 1) % 3] / .008) * std::sin(pi * p[(n + 2) % 3] / .008) *
			   std::cos(2 * pi * p[n] / .008);
	};
	compact.SetElectric(field);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int x = 0; x <= 4 * ratio[0]; ++x)
		{
			for (unsigned int y = 0; y <= 4 * ratio[1]; ++y)
			{
				for (unsigned int z = 0; z <= 4 * ratio[2]; ++z)
				{
					unsigned int pos[3] = {x, y, z};
					if (pos[n] == 4 * ratio[n])
					{
						continue;
					}
					double point[3];
					op->GetYeeCoords(n, pos, point, false);
					engine->SetVolt(n, pos,
									field(n, {point[0] * .001, point[1] * .001, point[2] * .001}) *
										op->GetEdgeLength(n, pos));
				}
			}
		}
	}
	for (unsigned int i = 0; i < 200; ++i)
	{
		compact.Step(dt);
		engine->IterateTS(1);
	}
	double peakE = 0, errorE = 0, peakH = 0, errorH = 0;
	for (unsigned int x = 1; x < 4 * ratio[0]; ++x)
	{
		for (unsigned int y = 1; y < 4 * ratio[1]; ++y)
		{
			for (unsigned int z = 1; z < 4 * ratio[2]; ++z)
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

static void CheckEnergy(unsigned int steps, Grid::Index ratio = {2, 2, 2})
{
	Grid grid({4, 4, 4}, {.002, .002, .002}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu, false, ratio);
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
	std::cout << ratio[0] << 'x' << ratio[1] << 'x' << ratio[2] << ": " << steps << " steps: relative invariant drift "
			  << error << '\n';
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

static double CavityFrequency(unsigned int cells, unsigned int component, Grid::Index ratio = {2, 2, 2})
{
	constexpr double length = .012;
	const double pi = std::acos(-1.);
	Grid grid({cells, cells, cells}, {length / cells, length / cells, length / cells},
			  {{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, epsilon, mu, false,
			  ratio);
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
	std::cout << "ratio " << ratio[0] << 'x' << ratio[1] << 'x' << ratio[2] << ": cavity base " << cells << ", E"
			  << component << ": frequency " << frequency << ", relative analytic error " << frequency / exact - 1
			  << '\n';
	return frequency / exact - 1;
}

static double PulseError(unsigned int cells, unsigned int propagation, Grid::Index ratio = {2, 2, 2})
{
	constexpr double length = .012;
	const Grid::Index dimensions{cells, cells, cells};
	const Grid::Point spacing{length / cells, length / cells, length / cells};
	const Grid::Box box{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}};
	Grid grid(dimensions, spacing, box, epsilon, mu, false, ratio);
	Grid reference(dimensions, spacing, box, epsilon, mu, true, ratio);
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
	std::cout << "ratio " << ratio[0] << 'x' << ratio[1] << 'x' << ratio[2] << ": pulse axis " << propagation
			  << ", base " << cells << ": peak relative error " << error / peak << '\n';
	Require(peak > 1e-3, "no pulse at reference receiver");
	return error / peak;
}

static void CheckAdapterRatio()
{
	ContinuousStructure csx;
	csx.GetGrid()->SetDeltaUnit(.001);
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 8; ++i)
		{
			csx.GetGrid()->AddDiscLine(n, i);
		}
	}
	Excitation excitation;
	excitation.SetupDiracPulse(30e9);
	std::unique_ptr<Operator_CartesianSubgrid> op(Operator_CartesianSubgrid::New({2, 2, 2, 6, 6, 6}, {1, 2, 8}));
	op->SetEngineType(Engine::BASIC);
	op->SetExcitationSignal(&excitation);
	Require(op->SetGeometryCSX(&csx), "adapter ratio geometry failed");
	int boundaries[6] = {0, 0, 0, 0, 0, 0};
	op->SetBoundaryCondition(boundaries);
	Require(op->CalcECOperator() == 0, "adapter ratio operator failed");
	Require(op->Grid().Ratio() == Grid::Index{1, 2, 8} && op->GetNumberCells() == 1472,
			"adapter dropped selected axis ratios");
	std::unique_ptr<Operator_CartesianSubgrid> multi(
		Operator_CartesianSubgrid::NewRegions({{1, 1, 1, 3, 4, 4}, {5, 1, 1, 7, 4, 4}}, {{2, 2, 4}, {4, 1, 2}}));
	multi->SetEngineType(Engine::BASIC);
	multi->SetExcitationSignal(&excitation);
	Require(multi->SetGeometryCSX(&csx), "multiple-region adapter geometry failed");
	multi->SetBoundaryCondition(boundaries);
	Require(multi->CalcECOperator() == 0 && multi->GetNumberCells() == 908 && multi->Grid().Regions().size() == 2,
			"native adapter dropped a refinement region");
}

static void CheckRatioGuards()
{
	for (const Grid::Index ratio : {Grid::Index{0, 2, 2}, {3, 2, 2}, {2, 32, 2}, {2, 2, ~0u}})
	{
		bool rejected = false;
		try
		{
			Grid bad({4, 4, 4}, {.001, .001, .001}, {{1, 1, 1}, {3, 3, 3}}, epsilon, mu, false, ratio);
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected, "invalid refinement ratio accepted");
	}
}

static std::vector<Grid::Region> SeparateRegions()
{
	return {{{{1, 2, 2}, {4, 6, 6}}, {2, 2, 4}}, {{{7, 2, 2}, {10, 6, 6}}, {4, 1, 2}}};
}

static void CheckMultipleRegions(unsigned int steps)
{
	const Grid::Index cells{12, 8, 8};
	const Grid::Point spacing{.001, .001, .001};
	auto regions = SeparateRegions();
	Grid grid(cells, spacing, epsilon, mu, regions);
	Require(grid.CellCount() == 1824 && grid.LatticeRatio() == Grid::Index{4, 2, 4},
			"multiple region count/lattice failed");
	double volume = 0;
	for (const auto &leaf : grid.Leaves())
	{
		const auto lattice = grid.LatticeRatio();
		Grid::Index root;
		for (unsigned int n = 0; n < 3; ++n)
		{
			root[n] = unsigned(leaf.lo[n]) / lattice[n];
		}
		Grid::Index expected{1, 1, 1};
		for (const auto &region : regions)
		{
			bool inside = true;
			for (unsigned int n = 0; n < 3; ++n)
			{
				inside = inside && root[n] >= region.box.start[n] && root[n] < region.box.stop[n];
			}
			if (inside)
			{
				expected = region.ratio;
			}
		}
		double local = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			Require(leaf.size[n] == int(lattice[n] / expected[n]), "a region inherited another region's refinement");
			local *= spacing[n] * leaf.size[n] / lattice[n];
		}
		volume += local;
	}
	Require(std::abs(volume / (.012 * .008 * .008) - 1) < 1e-12, "regions do not tile the physical domain");
	grid.SetElectric(RandomField);
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return RandomField(n, p) / 376.730313668; });
	double jump = 0;
	for (const auto &region : regions)
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
					for (unsigned int sample = 0; sample < 8; ++sample)
					{
						Grid::Point a;
						for (unsigned int n = 0; n < 3; ++n)
						{
							a[n] = spacing[n] * (region.box.start[n] +
												 (.17 + .07 * sample) * (region.box.stop[n] - region.box.start[n]));
						}
						a[normal] = spacing[normal] * (side ? region.box.stop[normal] : region.box.start[normal]);
						Grid::Point b = a;
						a[normal] -= 1e-12;
						b[normal] += 1e-12;
						jump = std::max(jump, std::abs(grid.Electric(tangent, a) - grid.Electric(tangent, b)));
					}
				}
			}
		}
	}
	const double residual = grid.CurlDivergenceResidual(RandomField);
	Require(jump < 1e-7 && residual < 1e-13, "multi-region trace or divergence of curl failed");
	// Region insertion order must not change propagation or interface ownership.
	std::reverse(regions.begin(), regions.end());
	Grid reversed(cells, spacing, epsilon, mu, regions);
	reversed.SetElectric(RandomField);
	reversed.SetMagnetic([](unsigned int n, const Grid::Point &p) { return RandomField(n, p) / 376.730313668; });
	const double dt = .8 * grid.StableTimestep();
	const double initial = grid.Energy(dt);
	double drift = 0;
	for (unsigned int step = 0; step < steps; ++step)
	{
		grid.Step(dt);
		if (step < 2000)
		{
			reversed.Step(dt);
			for (unsigned int n = 0; n < 3; ++n)
			{
				const Grid::Point p{.00817, .00327, .00461};
				Require(std::abs(grid.Electric(n, p) - reversed.Electric(n, p)) < 1e-12 &&
							std::abs(grid.Magnetic(n, p) - reversed.Magnetic(n, p)) < 1e-12,
						"region order changes sampled propagation");
			}
		}
		if (step % 100 == 0)
		{
			drift = std::max(drift, std::abs(grid.Energy(dt) / initial - 1));
			Require(std::isfinite(drift) && drift < 1e-8, "multi-region invariant is not conserved");
		}
	}
	std::cout << "multiple regions: cells " << grid.CellCount() << ", trace jump " << jump << ", D curl " << residual
			  << ", steps " << steps << ", invariant drift " << drift << '\n';
}

static void CheckManyRegionsAndBudget()
{
	std::vector<Grid::Region> regions;
	for (unsigned int x : {1u, 3u, 5u, 7u})
	{
		for (unsigned int y : {1u, 3u, 5u, 7u})
		{
			for (unsigned int z : {1u, 3u, 5u, 7u})
			{
				regions.push_back({{{x, y, z}, {x + 1, y + 1, z + 1}}, {2, 2, 2}});
			}
		}
	}
	Grid grid({10, 10, 10}, {.001, .001, .001}, epsilon, mu, regions);
	Require(grid.CellCount() == 1000 + 64 * 7, "64 independent regions have the wrong leaf count");
	Require(grid.CurlDivergenceResidual(RandomField) < 1e-13, "64-region divergence of curl failed");
	bool rejected = false;
	try
	{
		// Either block alone fits the guard. Their combined leaves do not.
		const std::vector<Grid::Region> oversized{{{{1, 1, 1}, {30, 31, 31}}, {2, 2, 4}},
												  {{{33, 1, 1}, {62, 31, 31}}, {2, 4, 4}}};
		Grid bad({64, 32, 32}, {.001, .001, .001}, epsilon, mu, oversized);
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "combined region allocation guard was not enforced");
}

static void CheckRegionGuards()
{
	const Grid::Region first{{{1, 2, 2}, {4, 6, 6}}, {2, 2, 4}};
	const std::vector<Grid::Region> invalid{
		{{{3, 2, 2}, {5, 6, 6}}, {2, 2, 2}},  // volume overlap
		{{{2, 3, 3}, {3, 4, 4}}, {2, 2, 2}},  // nested
		{{{4, 2, 2}, {6, 6, 6}}, {2, 2, 2}},  // face contact
		{{{4, 6, 6}, {6, 7, 7}}, {2, 2, 2}},  // corner contact
		{{{7, 2, 2}, {10, 6, 6}}, {3, 2, 2}}, // invalid ratio
		{{{0, 2, 2}, {1, 6, 6}}, {2, 2, 2}}	  // outer boundary
	};
	for (const auto &second : invalid)
	{
		bool rejected = false;
		try
		{
			Grid bad({12, 8, 8}, {.001, .001, .001}, epsilon, mu, {first, second});
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected, "invalid region interaction accepted");
	}
	for (const std::size_t count : {0u, 65u})
	{
		bool rejected = false;
		try
		{
			Grid bad({12, 8, 8}, {.001, .001, .001}, epsilon, mu, std::vector<Grid::Region>(count, first));
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected, "invalid region count accepted");
	}
}

static double MultipleRegionPulse(unsigned int cells, unsigned int propagation)
{
	constexpr double length = .012;
	const Grid::Point spacing{length / cells, length / cells, length / cells};
	Grid::Index ratioA{2, 2, 2}, ratioB{2, 2, 2};
	ratioA[propagation] = 2;
	ratioB[propagation] = 4;
	ratioB[(propagation + 2) % 3] = 1;
	Grid::Box a{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}};
	Grid::Box b = a;
	a.start[propagation] = cells / 6;
	a.stop[propagation] = cells / 3;
	b.start[propagation] = 2 * cells / 3;
	b.stop[propagation] = 5 * cells / 6;
	const std::vector<Grid::Region> regions{{a, ratioA}, {b, ratioB}};
	Grid grid({cells, cells, cells}, spacing, epsilon, mu, regions);
	Grid reference({cells, cells, cells}, spacing, epsilon, mu, regions, true);
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
	Grid::Point probe{.006, .006, .006};
	probe[propagation] = .011;
	const unsigned int steps = unsigned(1.1 * length * std::sqrt(epsilon * mu) / dt);
	double error = 0, peak = 0;
	for (unsigned int step = 0; step < steps; ++step)
	{
		grid.Step(dt);
		reference.Step(dt);
		const double value = reference.Electric(polarization, probe);
		peak = std::max(peak, std::abs(value));
		error = std::max(error, std::abs(value - grid.Electric(polarization, probe)));
	}
	Require(peak > 1e-3, "no packet beyond both refinement regions");
	std::cout << "multiple-region packet: axis " << propagation << ", base " << cells << ", error " << error / peak
			  << '\n';
	return error / peak;
}

int main(int argc, char **argv)
{
	try
	{
		if ((argc == 2 || argc == 5) && std::string(argv[1]) == "--mesh")
		{
			Grid::Index ratio{2, 2, 2};
			if (argc == 5)
			{
				for (unsigned int n = 0; n < 3; ++n)
				{
					const std::string token(argv[2 + n]);
					Require(token == "1" || token == "2" || token == "4" || token == "8" || token == "16",
							"invalid mesh export ratio");
					ratio[n] = std::stoul(token);
				}
			}
			Grid grid({8, 8, 8}, {.001, .001, .001}, {{2, 2, 2}, {6, 6, 6}}, epsilon, mu, false, ratio);
			std::cout << "{\"baseCells\":[8,8,8],\"spacingMm\":[1,1,1],\"refinedBoxMm\":[2,2,2,6,6,6],\"ratio\":["
					  << ratio[0] << ',' << ratio[1] << ',' << ratio[2] << "],\"leaves\":[";
			bool first = true;
			for (const auto &cell : grid.Leaves())
			{
				if (!first)
				{
					std::cout << ',';
				}
				first = false;
				std::cout << '[' << double(cell.lo[0]) / ratio[0] << ',' << double(cell.lo[1]) / ratio[1] << ','
						  << double(cell.lo[2]) / ratio[2] << ',' << double(cell.size[0]) / ratio[0];
				if (argc == 5)
				{
					std::cout << ',' << double(cell.size[1]) / ratio[1] << ',' << double(cell.size[2]) / ratio[2];
				}
				std::cout << ']';
			}
			std::cout << "],\"leafCount\":" << grid.CellCount() << ",\"accountedBytes\":" << grid.StorageBytes()
					  << "}\n";
			return 0;
		}
		CheckTraceAndStorage();
		const bool longRun = argc == 2 && std::string(argv[1]) == "--long";
		CheckMultipleRegions(longRun ? 100000 : 2000);
		CheckRegionGuards();
		CheckManyRegionsAndBudget();
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double a = MultipleRegionPulse(6, n);
			const double b = MultipleRegionPulse(12, n);
			Require(b < a, "multi-region packet does not converge");
		}
		const double multiFine = MultipleRegionPulse(24, 0);
		Require(multiFine < .1, "multi-region packet at resolution 24 differs by more than ten percent");
		const double multiResolved = MultipleRegionPulse(36, 0);
		Require(multiResolved < multiFine && multiResolved < .05,
				"resolved multi-region packet differs by more than five percent");
		const double anisotropicCoarse = CavityFrequency(6, 2, {2, 4, 8});
		const double anisotropicFine = CavityFrequency(12, 2, {2, 4, 8});
		Require(std::abs(anisotropicFine) < .005 && std::abs(anisotropicFine) < .75 * std::abs(anisotropicCoarse),
				"anisotropic cavity is not converging");
		if (longRun)
		{
			CheckEnergy(100000, {2, 4, 8});
		}
		CheckAdapterRatio();
		CheckRatioGuards();
		for (const Grid::Index ratio :
			 {Grid::Index{1, 1, 1}, {1, 2, 8}, {8, 1, 2}, {2, 8, 1}, {2, 4, 8}, {4, 4, 4}, {1, 1, 16}})
		{
			CheckTraceAndStorage(ratio);
			CheckEnergy(2000, ratio);
		}
		CheckUniformNative({1, 2, 4});
		CheckUniformNative({4, 1, 2});
		CheckUniformNative({2, 4, 1});
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Grid::Index ratio{1, 1, 1};
			ratio[axis] = 1;
			ratio[(axis + 1) % 3] = 2;
			ratio[(axis + 2) % 3] = 4;
			const double a = PulseError(6, axis, ratio);
			const double b = PulseError(12, axis, ratio);
			Require(b < .8 * a, "anisotropic packet fails to converge");
		}
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

/* SPDX-License-Identifier: GPL-3.0-or-later */
// Local-only spatial refinement against an independently qualified analytical reference.
#include "FDTD/cartesian_subgrid.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using Grid = CartesianSubgrid;
constexpr double epsilon = 8.85418781762e-12;
constexpr double mu = 1.2566370614359173e-6;
constexpr double unit = .01;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

struct Packet
{
	double angle, width, wavelength;
	bool magnetic;

	// Longitudinal packet centre -0.7, transverse ray offset +0.07.
	// Both TE E_z and TM H_z are divergence-free initial fields in XY.
	std::array<double, 3> ValueGradient(const Grid::Point &p) const
	{
		const double c = std::cos(angle), s = std::sin(angle);
		const double longitudinal = (c * p[0] + s * p[1]) / unit + .7;
		const double transverse = (-s * p[0] + c * p[1]) / unit - .07;
		const double wave = 2 * std::acos(-1.) / wavelength;
		const double envelope = std::exp(-std::pow(longitudinal / width, 2) -
			std::pow(transverse / .35, 2));
		const double value = envelope * std::cos(wave * longitudinal);
		const double a = -2 * longitudinal / (width * width) * value -
			wave * envelope * std::sin(wave * longitudinal);
		const double b = -2 * transverse / (.35 * .35) * value;
		return {value, (c * a - s * b) / unit, (s * a + c * b) / unit};
	}

	void Seed(Grid &grid, double timestep) const
	{
		grid.SetElectric([&](unsigned int n, const Grid::Point &p)
		{
			return !magnetic && n == 2 ? ValueGradient(p)[0] : 0.;
		});
		grid.SetMagnetic([&](unsigned int n, const Grid::Point &p)
		{
			const auto values = ValueGradient(p);
			if (magnetic)
			{
				return n == 2 ? values[0] / std::sqrt(mu / epsilon) : 0.;
			}
			// E is at t=0; H is at +dt/2 because this core updates E first.
			// First-half-step Faraday initialization, with O(dt^2) residual.
			if (n == 0)
			{
				return -.5 * timestep * values[2] / mu;
			}
			return n == 1 ? .5 * timestep * values[1] / mu : 0.;
		});
	}
};

static Grid MakeGrid(unsigned int cells, Grid::Index ratio, bool uniform, bool magnetic)
{
	Grid::Axes axes;
	for (auto &axis : axes)
	{
		for (unsigned int i = 0; i <= cells; ++i)
		{
			axis.push_back(unit * (-1.5 + 3. * i / cells));
		}
	}
	Grid::BoundaryConditions boundary{};
	if (magnetic)
	{
		// PMC caps admit the TM H_z mode without a spurious normal-H wall.
		boundary[4] = boundary[5] = Grid::BoundaryType::PMC;
	}
	return Grid(axes, epsilon, mu,
		{{{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, ratio}},
		uniform, {}, {}, boundary, {}, Grid::CurlLayout::Patterns, {}, 1200000);
}

// Integrate the existing local reconstruction over exactly the same physical
// cube in both grids. Its component basis is multilinear within each leaf,
// so the overlap midpoint integrates it exactly, without changing any field.
static Grid::Integral Average(const Grid &grid, const Grid::Point &centre,
	double halfWidth, bool magnetic)
{
	std::map<std::size_t, double> terms;
	double volume = 0;
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const auto lo = grid.LeafStart(leaf), lengths = grid.LeafLengths(leaf);
		Grid::Point midpoint;
		double overlap = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			const double a = std::max(lo[n], centre[n] - halfWidth);
			const double b = std::min(lo[n] + lengths[n], centre[n] + halfWidth);
			overlap *= std::max(0., b - a);
			midpoint[n] = .5 * (a + b);
		}
		if (overlap == 0)
		{
			continue;
		}
		const auto sample = magnetic ? grid.MagneticSample(2, midpoint) : grid.ElectricSample(2, midpoint);
		for (const auto &term : sample.weights)
		{
			terms[term.first] += overlap * term.second;
		}
		volume += overlap;
	}
	const double expected = std::pow(2 * halfWidth, 3);
	Require(std::abs(volume / expected - 1) < 1e-11, "averaging cube is not completely covered once");
	Grid::Integral result;
	result.magnetic = magnetic;
	for (const auto &term : terms)
	{
		result.weights.emplace_back(term.first, term.second / expected);
	}
	return result;
}

static void CheckSampling()
{
	const Packet packet{30 * std::acos(-1.) / 180, .25, 1.5, true};
	const Grid::Point point{unit * (std::cos(packet.angle) * .65 - std::sin(packet.angle) * .07),
		unit * (std::sin(packet.angle) * .65 + std::cos(packet.angle) * .07), .013 * unit};
	for (unsigned int n : {12u, 24u, 36u, 48u})
	{
		for (bool uniform : {false, true})
		{
			auto grid = MakeGrid(n, {2, 2, 2}, uniform, true);
			grid.SetMagnetic([](unsigned int component, const Grid::Point &p)
			{
				return component == 2 ? 1 + 2 * p[0] / unit - 3 * p[1] / unit : 0.;
			});
			const auto leaf = grid.LeafAt(point);
			const auto start = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
			const double expected = 1 + 2 * (start[0] + .5 * length[0]) / unit -
				3 * (start[1] + .5 * length[1]) / unit;
			Require(std::abs(grid.Magnetic(2, point) - expected) < 1e-12,
				"H point does not follow its documented transverse leaf-centre reconstruction");
			// Bounds +/-0.25 lie on every selected grid. The affine cube mean is1.
			const auto affine = Average(grid, {0, 0, 0}, .25 * unit, true);
			Require(std::abs(grid.Evaluate(affine) - 1) < 1e-12, "affine physical-volume average is wrong");
			grid.SetMagnetic([](unsigned int component, const Grid::Point &)
			{
				return component == 2 ? 1. : 0.;
			});
			const auto arbitrary = Average(grid, point, .08 * unit, true);
			Require(std::abs(grid.Evaluate(arbitrary) - 1) < 1e-12, "constant arbitrary-volume average is wrong");
			grid.SetElectric([](unsigned int component, const Grid::Point &p)
			{
				return component == 2 ? 1 + 2 * p[0] / unit - 3 * p[1] / unit : 0.;
			});
			Require(std::abs(grid.Evaluate(Average(grid, {0, 0, 0}, .25 * unit, false)) - 1) < 1e-12,
				"affine electric physical-volume average is wrong");
			std::cout << "SAMPLING_CONTROL base=" << n << " uniform=" << uniform
				<< " point_affine=" << expected << " analytic_point=" <<
				1 + 2 * point[0] / unit - 3 * point[1] / unit << '\n';
		}
	}
	std::cout << "PASS: unchanged H reconstruction and exact affine/constant volume controls" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		if (argc == 2 && std::string(argv[1]) == "--sampling-controls")
		{
			CheckSampling();
			return 0;
		}
		Require(argc == 12,
			"usage: local_packet N angle TE|TM rx ry rz width wavelength step_fraction csv baseline_dt");
		const unsigned int cells = std::stoul(argv[1]);
		Require(cells == 48 || cells == 72 || cells == 96, "unplanned local resolution");
		const double degrees = std::stod(argv[2]);
		Require(degrees == 30 || degrees == 45, "unplanned packet angle");
		const std::string polarization = argv[3];
		Require(polarization == "TE" || polarization == "TM", "expected TE or TM");
		Grid::Index ratio{static_cast<unsigned int>(std::stoul(argv[4])),
			static_cast<unsigned int>(std::stoul(argv[5])), static_cast<unsigned int>(std::stoul(argv[6]))};
		Require(ratio == Grid::Index{2, 2, 2} || ratio == Grid::Index{2, 4, 1}, "unplanned ratio");
		const Packet packet{degrees * std::acos(-1.) / 180, std::stod(argv[7]),
			std::stod(argv[8]), polarization == "TM"};
		const double fraction = std::stod(argv[9]);
		Require(packet.width >= .18 && packet.width <= .35 && packet.wavelength == 1.5 &&
			(fraction == .7 || fraction == .35), "unplanned source/time parameters");
		const double baselineStep = std::stod(argv[11]);
		Require(baselineStep > 0 && std::isfinite(baselineStep), "invalid frozen baseline timestep");
		auto local = MakeGrid(cells, ratio, false, packet.magnetic);
		// Reuse the accepted R03 step at N48 exactly. Finer trials scale that
		// common step with h, retaining a separate strict native CFL check.
		const double dt = cells == 48 ? baselineStep : baselineStep * (48. / cells);
		Require(dt <= local.StableTimestep(), "scaled timestep exceeds calculated local bound");
		packet.Seed(local, dt);
		const double initial = local.Energy(dt);
		Require(initial > 0, "initial packet has no energy");
		const double speed = 1 / std::sqrt(epsilon * mu);
		const unsigned int steps = static_cast<unsigned int>(std::ceil(2.4 * unit / (speed * dt)));
		std::vector<Grid::Point> points;
		for (double longitudinal : {-.6, 0., .65})
		{
			points.push_back({unit * (std::cos(packet.angle) * longitudinal - std::sin(packet.angle) * .07),
				unit * (std::sin(packet.angle) * longitudinal + std::cos(packet.angle) * .07), .013 * unit});
		}
		std::vector<Grid::Integral> averages;
		for (const auto &point : points)
		{
			averages.push_back(Average(local, point, .08 * unit, packet.magnetic));
		}
		std::ofstream output(argv[10]);
		Require(bool(output), "cannot create local history");
		output << std::setprecision(17) << "epoch,time_e,time_h";
		for (unsigned int p = 0; p < points.size(); ++p)
		{
			output << ",local_" << p;
		}
		for (unsigned int p = 0; p < points.size(); ++p)
		{
			output << ",avg_local_" << p;
		}
		output << '\n';
		double peak = 0, drift = 0;
		for (unsigned int step = 0; step < steps; ++step)
		{
			local.Step(dt);
			output << step + 1 << ',' << (step + 1) * dt << ',' << (step + 1.5) * dt;
			for (unsigned int p = 0; p < points.size(); ++p)
			{
				const double value = packet.magnetic ? local.Magnetic(2, points[p]) * std::sqrt(mu / epsilon) :
					local.Electric(2, points[p]);
				Require(std::isfinite(value), "nonfinite local packet field");
				output << ',' << value;
				if (p == 2)
				{
					peak = std::max(peak, std::abs(value));
				}
			}
			const double amplitude = packet.magnetic ? std::sqrt(mu / epsilon) : 1.;
			for (unsigned int p = 0; p < points.size(); ++p)
			{
				output << ',' << amplitude * local.Evaluate(averages[p]);
			}
			output << '\n';
			if (step % 16 == 0 || step + 1 == steps)
			{
				drift = std::max(drift, std::abs(local.Energy(dt) / initial - 1));
			}
		}
		output.close();
		Require(bool(output) && peak > .01 && drift < 1e-8, "local signal/invariant gate failed");
		std::cout << std::setprecision(17) << "LOCAL_PACKET_DIAGNOSTIC base=" << cells << " angle=" << degrees
			<< " polarization=" << polarization << " ratio=" << ratio[0] << ',' << ratio[1] << ',' << ratio[2]
			<< " leaves=" << local.CellCount() << " dt=" << dt << " baseline_dt=" << baselineStep
			<< " local_bound=" << local.StableTimestep() << " steps=" << steps << " peak=" << peak
			<< " local_drift=" << drift << " explicit_leaf_limit=1200000" << '\n';
		std::cout << "LOCAL_PACKET_FINISHED_ANALYTIC_REVIEW_REQUIRED" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << "LOCAL_PACKET_REJECTED: " << error.what() << std::endl;
		return 1;
	}
	return 0;
}

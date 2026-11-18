/* SPDX-License-Identifier: GPL-3.0-or-later */
// Device replay diagnostic only; CPU spatial point failures remain unchanged.
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_cartesiansubgrid.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
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
		uniform, {}, {}, boundary);
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
		Require(argc == 12, "usage: packet N angle TE|TM rx ry rz width wavelength step_fraction csv --reference|--metal|--cuda");
		const unsigned int cells = std::stoul(argv[1]);
		Require(cells == 12 || cells == 24 || cells == 36 || cells == 48, "unplanned base resolution");
		const double degrees = std::stod(argv[2]);
		Require(degrees >= 0 && degrees <= 90, "unsupported diagnostic angle");
		const std::string polarization = argv[3];
		Require(polarization == "TE" || polarization == "TM", "expected TE or TM");
		Grid::Index ratio{static_cast<unsigned int>(std::stoul(argv[4])),
			static_cast<unsigned int>(std::stoul(argv[5])), static_cast<unsigned int>(std::stoul(argv[6]))};
		Require(ratio == Grid::Index{1, 1, 1} || ratio == Grid::Index{2, 2, 2} ||
			ratio == Grid::Index{2, 4, 1}, "unplanned refinement ratio");
		const Packet packet{degrees * std::acos(-1.) / 180, std::stod(argv[7]),
			std::stod(argv[8]), polarization == "TM"};
		const double fraction = std::stod(argv[9]);
		Require(packet.width >= .18 && packet.width <= .35 && packet.wavelength == 1.5 &&
			(fraction == .7 || fraction == .35), "unplanned packet/time parameters");
		const std::string target = argv[11];
		Require(target == "--reference" || target == "--metal" || target == "--cuda",
			"expected explicit reference, Metal or CUDA target");
		if (target == "--cuda")
		{
			g_settings.SetGPUDevice(0);
		}
		const std::string expected = target == "--reference" ? "reference (" :
			target == "--metal" ? "Metal (" : "HIP (";
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target == "--reference" ? "reference" : "auto"));
		Require(bool(backend) && backend->GetName().find(expected) == 0,
			"requested device backend unavailable; no fallback allowed");
		std::cout << "PACKET_BACKEND name=" << backend->GetName() << std::endl;
		auto local = MakeGrid(cells, ratio, false, packet.magnetic);
		auto fine = MakeGrid(cells, ratio, true, packet.magnetic);
		const double dt = fraction * std::min(local.StableTimestep(), fine.StableTimestep());
		packet.Seed(local, dt);
		packet.Seed(fine, dt);
		auto localExecution = backend->CreateCartesianSubgrid(local);
		auto fineExecution = backend->CreateCartesianSubgrid(fine);
		Require(bool(localExecution) && bool(fineExecution), "selected device lacks subgrid execution");
		// Match the double CPU histories from R03. FP32 upload is part of the
		// measured backend difference, not silently folded into that reference.
		const double localInitial = local.Energy(dt), fineInitial = fine.Energy(dt);
		Require(localInitial > 0 && fineInitial > 0, "packet has no initial energy");
		const double speed = 1 / std::sqrt(epsilon * mu);
		const unsigned int steps = static_cast<unsigned int>(std::ceil(2.4 * unit / (speed * dt)));
		std::vector<Grid::Point> points;
		for (double longitudinal : {-.6, 0., .65})
		{
			points.push_back({unit * (std::cos(packet.angle) * longitudinal - std::sin(packet.angle) * .07),
				unit * (std::sin(packet.angle) * longitudinal + std::cos(packet.angle) * .07), .013 * unit});
		}
		std::ofstream output(argv[10]);
		Require(bool(output), "cannot create diagnostic histories");
		output << std::setprecision(17) << "epoch,time_e,time_h";
		for (unsigned int p = 0; p < points.size(); ++p)
		{
			output << ",local_" << p << ",fine_" << p;
		}
		std::vector<Grid::Integral> localAverage, fineAverage;
		for (unsigned int p = 0; p < points.size(); ++p)
		{
			localAverage.push_back(Average(local, points[p], .08 * unit, packet.magnetic));
			fineAverage.push_back(Average(fine, points[p], .08 * unit, packet.magnetic));
			output << ",avg_local_" << p << ",avg_fine_" << p;
		}
		output << '\n';
		double peak = 0, maximum = 0, localDrift = 0, fineDrift = 0;
		for (unsigned int step = 0; step < steps; ++step)
		{
			localExecution->Step(dt, 0);
			fineExecution->Step(dt, 0);
			// Every full epoch is observed for matched diagnostic histories. This
			// deliberate readback pattern is not a performance benchmark.
			localExecution->ReadFields(local);
			fineExecution->ReadFields(fine);
			output << step + 1 << ',' << (step + 1) * dt << ',' << (step + 1.5) * dt;
			for (unsigned int p = 0; p < points.size(); ++p)
			{
				const double a = packet.magnetic ? local.Magnetic(2, points[p]) * std::sqrt(mu / epsilon) :
					local.Electric(2, points[p]);
				const double b = packet.magnetic ? fine.Magnetic(2, points[p]) * std::sqrt(mu / epsilon) :
					fine.Electric(2, points[p]);
				Require(std::isfinite(a) && std::isfinite(b), "packet history became nonfinite");
				output << ',' << a << ',' << b;
				if (p == 2)
				{
					peak = std::max(peak, std::abs(b));
					maximum = std::max(maximum, std::abs(a - b));
				}
			}
			const double amplitude = packet.magnetic ? std::sqrt(mu / epsilon) : 1.;
			for (unsigned int p = 0; p < points.size(); ++p)
			{
				output << ',' << amplitude * local.Evaluate(localAverage[p]) <<
					',' << amplitude * fine.Evaluate(fineAverage[p]);
			}
			output << '\n';
			if (step % 16 == 0 || step + 1 == steps)
			{
				localDrift = std::max(localDrift, std::abs(local.Energy(dt) / localInitial - 1));
				fineDrift = std::max(fineDrift, std::abs(fine.Energy(dt) / fineInitial - 1));
			}
		}
		output.close();
		Require(bool(output) && peak > .01 && localDrift < 5e-4 && fineDrift < 5e-4,
			"missing transmitted packet or broken conservative evolution");
		std::cout << std::setprecision(17) << "PACKET_DIAGNOSTIC base=" << cells << " angle=" << degrees
			<< " polarization=" << polarization << " ratio=" << ratio[0] << ',' << ratio[1] << ',' << ratio[2]
			<< " leaves=" << local.CellCount() << " fine_leaves=" << fine.CellCount() << " dt=" << dt
			<< " local_bound=" << local.StableTimestep() << " fine_bound=" << fine.StableTimestep()
			<< " steps=" << steps << " peak=" << peak << " relative_peak_error=" << maximum / peak
			<< " local_drift=" << localDrift << " fine_drift=" << fineDrift << '\n';
		std::cout << "PACKET_DEVICE_STORAGE local_bytes=" << localExecution->StorageBytes()
			<< " fine_bytes=" << fineExecution->StorageBytes() << std::endl;
		std::cout << "DEVICE_DIAGNOSTIC_FINISHED_CPU_COMPARISON_REQUIRED" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << "DIAGNOSTIC_REJECTED: " << error.what() << std::endl;
		return 1;
	}
	return 0;
}

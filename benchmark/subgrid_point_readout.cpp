/* SPDX-License-Identifier: GPL-3.0-or-later */
// Read-only sampling geometry audit; no field propagation or solver changes.
#include "FDTD/cartesian_subgrid.h"
#include "tools/useful.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

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

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc == 4, "usage: geometry N 222|241 fresh.csv");
		const unsigned int cells = std::stoul(argv[1]);
		Require(cells == 48 || cells == 72 || cells == 96, "unplanned resolution");
		const std::string token = argv[2];
		Require(token == "222" || token == "241", "unplanned refinement ratio");
		const Grid::Index ratio = token == "222" ? Grid::Index{2, 2, 2} : Grid::Index{2, 4, 1};
		const double degrees = token == "222" ? 30 : 45;
		Grid::Axes axes;
		for (auto &axis : axes)
		{
			for (unsigned int i = 0; i <= cells; ++i)
			{
				axis.push_back(unit * (-1.5 + 3. * i / cells));
			}
		}
		Grid::BoundaryConditions boundary{};
		boundary[4] = boundary[5] = Grid::BoundaryType::PMC;
		Grid grid(axes, epsilon, mu,
			{{{{cells / 3, cells / 3, cells / 3}, {2 * cells / 3, 2 * cells / 3, 2 * cells / 3}}, ratio}},
			false, {}, {}, boundary, {}, Grid::CurlLayout::Patterns, {}, 1200000);
		// An independent affine field exposes the actual transverse-centre
		// reconstruction. It is a sampling control, not a Maxwell solution.
		grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
		{
			return n == 2 ? 1 + 2 * p[0] / unit - 3 * p[1] / unit + 4 * p[2] / unit : 0.;
		});
		std::ifstream existing(argv[3]);
		Require(!existing.good(), "refuse to overwrite geometry evidence");
		std::ofstream out(argv[3]);
		Require(out.good(), "unable to create geometry CSV");
		out << std::setprecision(17)
			<< "cells,rx,ry,rz,receiver,requested_x,requested_y,requested_z,centre_x,centre_y,centre_z,"
			<< "lo_x,lo_y,lo_z,size_x,size_y,size_z,leaf,leaves,affine_actual,affine_expected,integral_actual\n";
		const double angle = degrees * std::acos(-1.) / 180;
		unsigned int receiver = 0;
		for (double longitudinal : {-.6, 0., .65})
		{
			const Grid::Point point{
				unit * (std::cos(angle) * longitudinal - std::sin(angle) * .07),
				unit * (std::sin(angle) * longitudinal + std::cos(angle) * .07), .013 * unit};
			const auto leaf = grid.LeafAt(point);
			const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
			Grid::Point centre;
			for (unsigned int axis = 0; axis < 3; ++axis)
			{
				centre[axis] = lo[axis] + .5 * length[axis];
			}
			const double expected = 1 + 2 * centre[0] / unit - 3 * centre[1] / unit + 4 * point[2] / unit;
			const double actual = grid.Magnetic(2, point);
			const auto sample = grid.MagneticSample(2, point);
			const double integral = grid.Evaluate(sample);
			Require(std::abs(actual - expected) < 1e-12 && std::abs(integral - expected) < 1e-12,
				"native H reconstruction does not match its documented affine geometry");
			out << cells << ',' << ratio[0] << ',' << ratio[1] << ',' << ratio[2] << ',' << receiver++;
			for (const auto &p : {point, centre, lo, length})
			{
				for (double value : p)
				{
					out << ',' << value;
				}
			}
			out << ',' << leaf << ',' << grid.CellCount() << ',' << actual << ',' << expected << ',' << integral << '\n';
		}
		out.close();
		Require(out.good(), "geometry CSV write failed");
		std::cout << "NATIVE_READOUT_GEOMETRY_AUDIT_FINISHED cells=" << cells
			<< " ratio=" << token << " leaves=" << grid.CellCount() << " affine_controls=3 no_fdtd=1" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

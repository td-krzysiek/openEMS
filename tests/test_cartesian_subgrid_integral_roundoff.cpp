/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using Grid = CartesianSubgrid;
static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}
static void Check(unsigned int normal, double offset)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			axes[n].push_back((n == normal ? offset : 0) + .001 * i);
		}
	}
	Grid::Index ratio{1, 1, 1};
	ratio[normal] = 2;
	Grid grid(axes, 8.85418781762e-12, 1.2566370614359173e-6,
		{{{{1, 1, 1}, {3, 3, 3}}, ratio}});
	const auto direction = (normal + 1) % 3, cross = (normal + 2) % 3;
	Grid::Point start{.002, .002, .002};
	start[normal] = axes[normal][2];
	start[direction] = .001;
	start[cross] = .002;
	auto stop = start;
	stop[direction] = .002;
	const auto exact = grid.VoltagePath(start, stop);
	Require(exact.weights.size() == 1, "exact local edge has spurious voltage support");
	for (int sign : {-1, 1})
	{
		auto a = start, b = stop;
		for (unsigned int ulp = 0; ulp < 2; ++ulp)
		{
			a[normal] = std::nextafter(a[normal], sign < 0 ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity());
		}
		b[normal] = a[normal];
		const auto rounded = grid.VoltagePath(a, b);
		Require(rounded.weights == exact.weights,
			"machine-roundoff boundary lookup creates extra voltage-path support");
	}
	// A genuinely interior location retains both basis contributions.
	auto interiorA = start, interiorB = stop;
	interiorA[normal] -= 1e-10;
	interiorB[normal] = interiorA[normal];
	const auto interior = grid.VoltagePath(interiorA, interiorB);
	Require(interior.weights.size() > 1, "roundoff repair snapped a resolved interior path");
	for (const auto &term : interior.weights)
	{
		Require(term.second > 0, "interior voltage path has a negative basis weight");
	}
	// At a refinement face the chosen positive-side magnetic chain needs all
	// fine cuts, even when a transverse coordinate rounds just below the face.
	Grid::Index tensorRatio{4, 4, 4};
	tensorRatio[normal] = 2;
	Grid tensor(axes, 8.85418781762e-12, 1.2566370614359173e-6,
		{{{{1, 1, 1}, {3, 3, 3}}, tensorRatio}});
	Grid::Point contourA{}, contourB{};
	for (unsigned int n = 0; n < 3; ++n)
	{
		contourA[n] = axes[n][1];
		contourB[n] = axes[n][n == normal ? 1 : 3];
	}
	const auto exactContour = tensor.CurrentContour(normal, contourA, contourB);
	Require(exactContour.weights.size() > 8, "contour control lacks refined magnetic support");
	for (int sign : {-1, 1})
	{
		auto a = contourA, b = contourB;
		for (unsigned int ulp = 0; ulp < 2; ++ulp)
		{
			a[normal] = std::nextafter(a[normal], sign < 0 ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity());
		}
		b[normal] = a[normal];
		Require(tensor.CurrentContour(normal, a, b).weights == exactContour.weights,
			"roundoff boundary lookup misses fine magnetic quadrature cuts");
	}
	std::cout << "INTEGRAL_ROUNDOFF normal=" << normal << " offset=" << offset << " interior=" << interior.weights.size() << std::endl;
}
int main()
{
	try
	{
		for (unsigned int normal = 0; normal < 3; ++normal)
		{
			for (double offset : {.034, -.038, 0.})
			{
				Check(normal, offset);
			}
		}
		std::cout << "PASS: boundary-roundoff integrals preserve exact support and resolved interior weights" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

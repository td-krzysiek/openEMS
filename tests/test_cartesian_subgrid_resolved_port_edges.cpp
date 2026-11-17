/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "tools/constants.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

using Grid = CartesianSubgrid;
static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}
static double RelativeTraceError(const Grid::Integral &expected, const Grid::Integral &actual)
{
	std::map<std::size_t, double> difference;
	double norm = 0, error = 0;
	for (const auto &term : expected.weights)
	{
		Require(std::isfinite(term.second), "nonfinite general voltage weight");
		difference[term.first] += term.second;
		norm += std::abs(term.second);
	}
	for (const auto &term : actual.weights)
	{
		Require(std::isfinite(term.second), "nonfinite resolved voltage weight");
		difference[term.first] -= term.second;
	}
	for (const auto &term : difference)
	{
		error += std::abs(term.second);
	}
	Require(norm > 0 && std::isfinite(norm) && std::isfinite(error), "invalid voltage trace norm");
	return error / norm;
}
static void Control(unsigned int axis, bool graded, bool refined, bool unequal, bool sheet)
{
	Grid::Axes base, fine;
	for (unsigned int n = 0; n < 3; ++n)
	{
		base[n] = {0, graded ? .0007 + .0001 * n : .001, graded ? .0018 + .0001 * n : .002,
			graded ? .0031 + .0001 * n : .003, .004};
		for (std::size_t k = 0; k < 4; ++k)
		{
			fine[n].push_back(base[n][k]);
			fine[n].push_back(base[n][k] + (unequal ? .37 : .5) * (base[n][k + 1] - base[n][k]));
		}
		fine[n].push_back(base[n].back());
	}
	const Grid::Region region = refined ? Grid::Region{{{1, 1, 1}, {3, 3, 3}}, {2, 2, 2}} :
		Grid::Region{{{2, 2, 2}, {6, 6, 6}}, {1, 1, 1}};
	Grid grid(refined ? base : fine, EPS0 * 3.556, MUE0, {region}, false, {}, {}, {}, {},
		Grid::CurlLayout::Patterns, refined && unequal ? fine : Grid::Axes{});
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return .2 + n + .5 * p[0];
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return std::cos(1000 * p[0] + n) / 376.730313668;
	});
	const auto before = grid.DeviceStateData();
	Grid::Point a, b;
	for (unsigned int n = 0; n < 3; ++n)
	{
		a[n] = fine[n][refined ? 3 : 1];
		b[n] = fine[n][refined ? 5 : 7];
	}
	if (sheet)
	{
		a[(axis + 2) % 3] = b[(axis + 2) % 3] = fine[(axis + 2) % 3][4];
	}
	const auto resolved = grid.PrepareDistributedResistor(a, b, axis, 50);
	const auto oldAPI = grid.DistributedResistor(a, b, axis, 50);
	const auto reverse = grid.PrepareDistributedResistor(b, a, axis, 50);
	Require(resolved.size() == oldAPI.size() && resolved.size() == reverse.size() && !resolved.empty(),
		"resolved descriptor count changed");
	double maximum = 0;
	std::size_t exact = 0;
	for (std::size_t i = 0; i < resolved.size(); ++i)
	{
		const auto &entry = resolved[i];
		Require(entry.resistor.start == oldAPI[i].start && entry.resistor.stop == oldAPI[i].stop &&
			entry.resistor.resistance == oldAPI[i].resistance && entry.edge == reverse[i].edge &&
			entry.resistor.start == reverse[i].resistor.start && entry.resistor.stop == reverse[i].resistor.stop &&
			entry.resistor.resistance == reverse[i].resistor.resistance &&
			entry.leaf == reverse[i].leaf && entry.first == reverse[i].first && entry.second == reverse[i].second,
			"prepared support changes normalization, endpoints or reversed input");
		const auto slow = grid.VoltagePath(entry.resistor.start, entry.resistor.stop);
		const auto direct = grid.LeafEdgeVoltage(entry.leaf, axis, entry.first, entry.second);
		Require(direct.weights.size() == 1 && direct.weights.front().first == entry.edge,
			"resolved independent edge points at another trace");
		const double error = RelativeTraceError(slow, direct);
		maximum = std::max(maximum, error);
		exact += slow.weights == direct.weights;
		// Keep the original2e-14 numerical bound, now on the complete functional
		// (a stronger L1 check). General quadrature can retain tiny extra terms
		// at rounded cut endpoints; requiring identical support was not valid.
		auto wrong = direct;
		wrong.weights.push_back({entry.edge ? 0 : 1, .01 * direct.weights.front().second});
		Require(RelativeTraceError(slow, wrong) > 2e-14, "wrong support escapes the trace-error guard");
		wrong = direct;
		wrong.weights.front().second *= 1.01;
		Require(RelativeTraceError(slow, wrong) > 2e-14, "wrong main weight escapes the trace-error guard");
	}
	// Coordinate subtraction in the old general quadrature and the stored
	// local length can differ by roundoff. State that explicitly; retain a
	// tight path bound as well as all existing field/power acceptance gates.
	Require(maximum < 2e-14, "direct local voltage differs from general integration");
	const auto after = grid.DeviceStateData();
	Require(before.electric == after.electric && before.magnetic == after.magnetic,
		"preparing port descriptors modified fields");
	std::cout << "RESOLVED_PORT axis=" << axis << " graded=" << graded << " refined=" << refined
		<< " unequal=" << unequal << " sheet=" << sheet << " edges=" << resolved.size()
		<< " byte_exact_weights=" << exact << " maximum_relative_weight_error=" << maximum
		<< " fields_unchanged=1" << std::endl;
}
int main()
{
	try
	{
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool graded : {false, true})
			{
				for (bool refined : {false, true})
				{
					for (bool unequal : {false, true})
					{
						for (bool sheet : {false, true})
						{
							Control(axis, graded, refined, unequal, sheet);
						}
					}
				}
			}
		}
		std::cout << "PASS: resolved local port support, normalization, roundoff bounds and unchanged fields" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>

using Grid = CartesianSubgrid;
static constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

static Grid Fixture(unsigned int mode, Grid::CurlLayout layout)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		axes[n].push_back(0);
		for (unsigned int i = 0; i < 10; ++i)
		{
			axes[n].push_back(axes[n].back() + .001 * (1 + .1 * n + .03 * i));
		}
	}
	std::vector<Grid::Region> regions;
	if (mode == 0)
	{
		regions.push_back({{{2, 2, 3}, {4, 8, 7}}, {1, 1, 1}});
	}
	if (mode != 0)
	{
		regions.push_back({{{2, 2, 3}, {4, 8, 7}}, {2, 4, 2}});
	}
	if (mode == 2)
	{
		regions.push_back({{{6, 3, 3}, {8, 7, 7}}, {4, 1, 2}});
	}
	const auto material = [&](const Grid::Point &, const Grid::Point &hi)
	{
		return Grid::Medium{epsilon * (hi[0] <= axes[0][5] ? 3.556 : 1), mu};
	};
	const auto pec = [&](unsigned int, const Grid::Point &lo, const Grid::Point &hi)
	{
		return mode == 4 && lo[2] == axes[2][1] && hi[2] == axes[2][1];
	};
	const auto opaque = [&](const Grid::Point &lo, const Grid::Point &hi)
	{
		return mode == 5 && lo[2] >= axes[2][1] && hi[2] <= axes[2][2];
	};
	return Grid(axes, epsilon, mu, regions, mode == 3, material, pec, {}, opaque, layout);
}

static void Control(unsigned int mode, Grid::CurlLayout layout)
{
	const auto grid = Fixture(mode, layout);
	const auto expected = grid.DeviceCurls(false);
	std::vector<std::pair<std::size_t, double>> row;
	row.reserve(32);
	std::size_t terms = 0;
	for (std::size_t face = 0; face < grid.FaceCount(); ++face)
	{
		grid.CurlRow(face, row);
		const auto pattern = expected.curl.pattern.at(face), base = expected.curl.base.at(face);
		const auto begin = expected.curl.start.at(pattern), end = expected.curl.start.at(pattern + 1);
		Require(row.size() == end - begin, "direct curl row changed its length");
		for (std::size_t k = 0; k < row.size(); ++k)
		{
			const auto &term = row[k];
			const float high = static_cast<float>(term.second);
			const float low = static_cast<float>(term.second - double(high));
			Require(term.first == base + expected.curl.offset[begin + k] &&
				std::memcmp(&high, &expected.curl.weight[begin + k], sizeof(float)) == 0 &&
				std::memcmp(&low, &expected.curl.residual[begin + k], sizeof(float)) == 0,
				"direct normalized row changed a canonical coefficient or order");
		}
		terms += row.size();
	}
	// A caller may reuse one buffer for arbitrary rows, including empty rows.
	for (std::size_t face = grid.FaceCount(); face > 0; --face)
	{
		grid.CurlRow(face - 1, row);
		const auto pattern = expected.curl.pattern.at(face - 1);
		Require(row.size() == expected.curl.start[pattern + 1] - expected.curl.start[pattern],
			"reused direct row retains a previous row's entries");
	}
	const auto retained = row;
	try
	{
		grid.CurlRow(grid.FaceCount(), row);
		throw std::runtime_error("invalid curl row was accepted");
	}
	catch (const std::out_of_range &)
	{
		Require(row == retained, "rejected row access changed caller storage");
	}
	std::cout << "DIRECT_CURL_ROWS mode=" << mode
		<< " spatial=" << (layout == Grid::CurlLayout::SpatialBlocks)
		<< " rows=" << grid.FaceCount() << " terms=" << terms << " coefficients_exact=1" << std::endl;
}

int main()
{
	try
	{
		for (auto layout : {Grid::CurlLayout::Patterns, Grid::CurlLayout::SpatialBlocks})
		{
			for (unsigned int mode = 0; mode < 6; ++mode)
			{
				Control(mode, layout);
			}
		}
		std::cout << "PASS: direct canonical curl rows, reusable storage and bounds" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

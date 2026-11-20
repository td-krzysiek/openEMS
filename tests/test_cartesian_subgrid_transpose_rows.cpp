/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_child_knots_fixture.h"
#include <cstring>
#include <iostream>
#include <stdexcept>

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

static Grid Make(unsigned int mode, Grid::CurlLayout layout)
{
	const auto base = Base();
	Grid::Index ratio = mode == 7 ? Grid::Index{16, 2, 1} : Grid::Index{2, 4, 2};
	if (mode == 9)
	{
		ratio = {2, 16, 1};
	}
	if (mode == 10)
	{
		ratio = {1, 2, 16};
	}
	if (mode == 11)
	{
		ratio = {4, 4, 4};
	}
	std::vector<Grid::Region> regions{{{{2, 2, 2}, {3, 4, 4}}, ratio}};
	if (mode == 0)
	{
		regions = {{{{2, 2, 2}, {3, 4, 4}}, {1, 1, 1}}};
		ratio = {1, 1, 1};
	}
	if (mode == 2 || mode == 6)
	{
		regions.push_back({{{5, 3, 3}, {6, 6, 6}}, {4, 1, 2}});
		ratio = {4, 4, 2};
	}
	const auto medium = [&](const Grid::Point &, const Grid::Point &hi)
	{
		return Grid::Medium{epsilon * (hi[0] <= base[0][4] ? 3.25 : 1), mu};
	};
	const auto pec = [&](unsigned int, const Grid::Point &lo, const Grid::Point &hi)
	{
		return mode == 4 && lo[2] == base[2][1] && hi[2] == base[2][1];
	};
	const auto opaque = [&](const Grid::Point &lo, const Grid::Point &hi)
	{
		return mode == 5 && lo[2] >= base[2][0] && hi[2] <= base[2][1];
	};
	Grid::BoundaryConditions boundaries{};
	if (mode == 8)
	{
		boundaries.fill(Grid::BoundaryType::PMC);
	}
	return Grid(base, epsilon, mu, regions, mode == 3, medium, pec, boundaries,
		opaque, layout, mode == 6 || mode == 7 || mode >= 9 ? Children(base, ratio) : Grid::Axes{});
}

static void Control(unsigned int mode, Grid::CurlLayout layout)
{
	auto grid = Make(mode, layout);
	// The old expanded C^T is built independently by scanning every C row.
	const auto oracle = grid.DeviceCurls(false);
	std::vector<std::pair<std::size_t, double>> row;
	row.reserve(32);
	std::size_t terms = 0;
	for (std::size_t id = 0; id < grid.EdgeCount(); ++id)
	{
		grid.CurlTransposeRow(id, row);
		const auto pattern = oracle.transpose.pattern.at(id);
		const auto base = oracle.transpose.base.at(id);
		const auto begin = oracle.transpose.start.at(pattern);
		const auto end = oracle.transpose.start.at(pattern + 1);
		if (row.size() != end - begin)
		{
			std::cerr << "mode=" << mode << " edge=" << id << " actual=" << row.size()
				<< " expected=" << end - begin << std::endl;
			throw std::runtime_error("implicit transpose changed support size");
		}
		for (std::size_t k = 0; k < row.size(); ++k)
		{
			const auto &term = row[k];
			const float high = static_cast<float>(term.second);
			const float low = static_cast<float>(term.second - double(high));
			Require(term.first == base + oracle.transpose.offset[begin + k] &&
				std::memcmp(&high, &oracle.transpose.weight[begin + k], sizeof(float)) == 0 &&
				std::memcmp(&low, &oracle.transpose.residual[begin + k], sizeof(float)) == 0,
				"implicit transpose changed canonical coefficient or order");
		}
		terms += row.size();
	}
	const auto retained = row;
	try
	{
		grid.CurlTransposeRow(grid.EdgeCount(), row);
		throw std::runtime_error("invalid transpose edge accepted");
	}
	catch (const std::out_of_range &)
	{
		Require(row == retained, "invalid transpose access changed caller buffer");
	}
	grid.UpdateElectric(.1 * grid.StableTimestep());
	try
	{
		grid.CurlTransposeRow(0, row);
		throw std::runtime_error("incomplete timestep accepted");
	}
	catch (const std::invalid_argument &)
	{
		Require(row == retained, "phase rejection changed caller buffer");
	}
	std::cout << "IMPLICIT_TRANSPOSE mode=" << mode
		<< " spatial=" << (layout == Grid::CurlLayout::SpatialBlocks)
		<< " rows=" << grid.EdgeCount() << " terms=" << terms
		<< " coefficients_exact=1" << std::endl;
}

int main()
{
	try
	{
		for (auto layout : {Grid::CurlLayout::Patterns, Grid::CurlLayout::SpatialBlocks})
		{
			for (unsigned int mode = 0; mode < 12; ++mode)
			{
				Control(mode, layout);
			}
		}
		std::cout << "PASS: implicit canonical transpose support and rejected access" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

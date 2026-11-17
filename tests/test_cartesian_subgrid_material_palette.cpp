/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/subgrid_material_palette.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;
using Table = CartesianSubgridMaterialPalette<Grid::Medium>;
static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}
template <class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid palette lifecycle/access accepted");
}
static bool Exact(const Grid::Medium &a, const Grid::Medium &b)
{
	return std::memcmp(&a.epsilon, &b.epsilon, sizeof(double)) == 0 &&
		std::memcmp(&a.permeability, &b.permeability, sizeof(double)) == 0;
}
static void TableControl(unsigned int mode)
{
	const std::size_t count = mode == 2 ? 5000 : 10000;
	Table table;
	table.Reserve(count);
	Reject([&] { table.Reserve(count); });
	Reject([&] { table.Finish(); });
	std::vector<Grid::Medium> oracle;
	for (std::size_t i = 0; i < count; ++i)
	{
		double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
		if (mode == 1)
		{
			epsilon = i % 3 == 0 ? epsilon : std::nextafter(epsilon, 1.);
			mu = i % 5 == 0 ? std::nextafter(mu, 0.) : mu;
		}
		else if (mode == 2)
		{
			epsilon *= 1 + .0001 * i;
		}
		oracle.push_back({epsilon, mu});
		table.Append(oracle.back());
		Require(Exact(table.at(i), oracle.back()), "append changed material bits");
	}
	table.Finish();
	Require(table.Size() == count && table.Dense() == (mode == 2), "wrong palette/dense mode");
	Require(table.UniqueCount() == (mode == 0 ? 1u : mode == 1 ? 4u : 0u), "wrong exact palette size");
	if (mode == 0)
	{
		Require(table.StorageBytes() == sizeof(Grid::Medium), "homogeneous table allocated per-leaf IDs");
	}
	else if (mode == 1)
	{
		Require(table.StorageBytes() <= count * sizeof(std::uint32_t) + 4 * sizeof(Grid::Medium),
			"palette retained dense media or setup dictionary");
	}
	else
	{
		Require(table.StorageBytes() == count * sizeof(Grid::Medium), "dense fallback retained palette IDs");
	}
	auto copy = table;
	auto moved = std::move(copy);
	for (std::size_t i = 0; i < count; ++i)
	{
		Require(Exact(table[i], oracle[i]) && Exact(moved[i], oracle[i]), "copy/move changed material ownership");
	}
	Reject([&] { table.Append(oracle.front()); });
	Reject([&] { table.Finish(); });
	Reject([&] { table.at(count); });
	std::cout << "MATERIAL_PALETTE mode=" << mode << " leaves=" << count
		<< " bytes=" << table.StorageBytes() << " exact_copy=1 lifecycle_guards=5\n";
}
static void GridControl()
{
	Grid::Axes lines;
	for (auto &axis : lines)
	{
		axis = {0, .001, .002, .003, .004, .005, .006, .007, .008};
	}
	const std::vector<Grid::Region> regions{{{{2, 2, 2}, {6, 6, 6}}, {2, 4, 2}}};
	for (bool varying : {false, true})
	{
		const auto material = [&](const Grid::Point &lo, const Grid::Point &)
		{
			return Grid::Medium{8.85418781762e-12 * (varying ? 1 + 19 * lo[0] + 29 * lo[1] + 43 * lo[2] : 3.556),
				1.2566370614359173e-6 * (varying && lo[2] >= .004 ? 1.7 : 1)};
		};
		Grid grid(lines, 8.85418781762e-12, 1.2566370614359173e-6, regions, false, material);
		for (std::size_t i = 0; i < grid.CellCount(); ++i)
		{
			const auto lo = grid.LeafStart(i), len = grid.LeafLengths(i);
			Grid::Point hi;
			for (unsigned int n = 0; n < 3; ++n)
			{
				hi[n] = lo[n] + len[n];
			}
			Require(Exact(grid.LeafMedium(i), material(lo, hi)), "grid material sample changed");
		}
		Reject([&] { grid.LeafMedium(grid.CellCount()); });
		std::cout << "MATERIAL_GRID varying=" << varying << " leaves=" << grid.CellCount() << " exact_bits=1\n";
	}
}
int main()
{
	try
	{
		for (unsigned int mode = 0; mode < 3; ++mode)
		{
			TableControl(mode);
		}
		GridControl();
		std::cout << "PASS: exact bounded material palette, dense fallback and grid samples\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

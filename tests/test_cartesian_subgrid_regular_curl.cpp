/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

using Grid = CartesianSubgrid;
static void Require(bool value, const char *message)
{
	if (!value)
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
	catch (const std::invalid_argument &)
	{
		return;
	}
	throw std::runtime_error("invalid structured metadata or missing metric view was admitted");
}
static void IndependentFaceOracle()
{
	using Regular = CartesianSubgridRegularCurl;
	Regular::Block block;
	block.cells = {1, 1, 1};
	block.length = {{{1}, {2}, {3}}};
	block.edge = {{{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}}};
	std::vector<double> mass(12, 2), reluctance(6, 3);
	const CartesianSubgridMetricView metrics{mass.data(), reluctance.data(), mass.size(), reluctance.size()};
	using Row = std::vector<std::pair<std::size_t, double>>;
	std::vector<Row> rows{
		{{4, 3}, {5, -3}, {8, -4.5}, {9, 4.5}},
		{{6, 3}, {7, -3}, {10, -4.5}, {11, 4.5}},
		{{0, -1.5}, {1, 1.5}, {8, 4.5}, {10, -4.5}},
		{{2, -1.5}, {3, 1.5}, {9, 4.5}, {11, -4.5}},
		{{0, 1.5}, {2, -1.5}, {4, -3}, {6, 3}},
		{{1, 1.5}, {3, -1.5}, {5, -3}, {7, 3}}};
	// Do not round a supplied coefficient into a convenient geometric pattern.
	rows[4][0].second = std::nextafter(rows[4][0].second, 2.);
	CartesianSubgridCurl curl;
	curl.ConfigureBlocks({block}, {{0, 0}}, 12);
	CartesianSubgridCurl::Builder builder(curl, 6, 12, metrics);
	for (const auto &row : rows)
	{
		builder.Append(row.begin(), row.end());
	}
	builder.Finish();
	Require(curl.RegularRowCount() == 5, "one-ULP exception was silently quantized");
	for (std::size_t row = 0; row < 6; ++row)
	{
		Row actual;
		curl.ForEach(row, [&](std::size_t edge, double value) { actual.push_back({edge, value}); }, metrics);
		Require(actual == rows[row], "independent face orientation or exceptional coefficient changed");
	}
	Reject([&] { curl.ForEach(0, [](std::size_t, double) {}); });
	std::vector<double> e(12), h(6), expected(6, 0), transposed(12, 0), actualTranspose(12);
	for (std::size_t i = 0; i < e.size(); ++i)
	{
		e[i] = std::sin(.73 * i);
	}
	for (std::size_t row = 0; row < 6; ++row)
	{
		for (const auto &term : rows[row])
		{
			expected[row] += term.second * e[term.first];
		}
	}
	curl.Apply(e, h, metrics);
	Require(h == expected, "structured gather differs from the independent face matrix");
	for (std::size_t row = 0; row < 6; ++row)
	{
		for (const auto &term : rows[row])
		{
			transposed[term.first] += term.second * h[row];
		}
	}
	curl.ApplyTranspose(h, actualTranspose, metrics);
	Require(transposed == actualTranspose, "structured transpose changed order or work pairing");
	Reject([&] { CartesianSubgridCurl::Builder invalid(curl, 6, 12); });
	Require(curl.RegularRowCount() == 5, "invalid rebuild destroyed prior structured rows");
	Row retained;
	curl.ForEach(4, [&](std::size_t edge, double value) { retained.push_back({edge, value}); }, metrics);
	Require(retained == rows[4], "invalid rebuild destroyed prior exceptional coefficients");
	Regular invalid;
	block.edge[0][0] = 12;
	Reject([&] { invalid.Configure({block}, {{0, 0}}, 12); });
}
static void Compare(const std::string &name)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= (name == "many" ? 10u : 8u); ++i)
		{
			axes[n].push_back(.001 * i * (1 + .1 * n) + (name == "graded" ? 1e-6 * i * i : 0));
		}
	}
	std::vector<Grid::Region> regions{{{{2, 2, 2}, {6, 6, 6}}, name == "identity" ? Grid::Index{1, 1, 1} : Grid::Index{2, 2, 2}}};
	if (name == "multiple")
	{
		regions = {{{{1, 1, 1}, {3, 3, 3}}, {2, 4, 1}}, {{{5, 5, 5}, {7, 7, 7}}, {4, 1, 2}}};
	}
	if (name == "many")
	{
		regions.clear();
		for (unsigned int x : {1, 3, 5, 7})
		{
			for (unsigned int y : {1, 3, 5, 7})
			{
				for (unsigned int z : {1, 3, 5, 7})
				{
					regions.push_back({{{x, y, z}, {x + 1, y + 1, z + 1}}, {2, 2, 2}});
				}
			}
		}
	}
	if (name == "mixed16")
	{
		regions = {{{{2, 2, 2}, {3, 3, 3}}, {16, 8, 2}}};
	}
	Grid::MediumSampler material;
	if (name == "material")
	{
		material = [&](const Grid::Point &lo, const Grid::Point &)
		{
			return Grid::Medium{8.85418781762e-12 * (lo[0] >= axes[0][4] ? 2.3 : 1),
				1.2566370614359173e-6};
		};
	}
	auto make = [&](Grid::CurlLayout layout)
	{
		return Grid(axes, 8.85418781762e-12, 1.2566370614359173e-6, regions,
			name == "uniform-fine", material, {}, {}, {}, layout);
	};
	Reject([&] { make(static_cast<Grid::CurlLayout>(99)); });
	auto reference = make(Grid::CurlLayout::Patterns), blocks = make(Grid::CurlLayout::SpatialBlocks);
	Require(blocks.StableTimestep() == reference.StableTimestep(), "block storage changed the coupled timestep");
	const auto a = reference.DeviceData(), b = blocks.DeviceData();
	Require(a.curlStart == b.curlStart && a.curlColumn == b.curlColumn && a.curlWeight == b.curlWeight && a.curlLow == b.curlLow &&
		a.transposeStart == b.transposeStart && a.transposeColumn == b.transposeColumn &&
		a.transposeWeight == b.transposeWeight && a.transposeLow == b.transposeLow &&
		a.blockStart == b.blockStart && a.blockEdge == b.blockEdge && a.lower == b.lower,
		"spatial blocks changed a device curl, transpose or interface factor");
	Require(blocks.StructuredCurlStorage().regularRows > blocks.FaceCount() / 4, "block path was silently bypassed");
	Require(reference.StructuredCurlStorage().regularRows == 0, "CPU reference changed its representation");
	for (Grid *grid : {&reference, &blocks})
	{
		grid->SetElectric([](unsigned int n, const Grid::Point &p) { return std::sin(731 * p[(n + 1) % 3]); });
		grid->SetMagnetic([](unsigned int n, const Grid::Point &p) { return .001 * std::cos(317 * p[n]); });
		grid->SetSources({{0, .25}});
	}
	const double dt = .5 * reference.StableTimestep();
	for (unsigned int i = 0; i < 128; ++i)
	{
		reference.Step(dt, .2 * std::sin(.1 * i));
		blocks.Step(dt, .2 * std::sin(.1 * i));
		Require(reference.Energy(dt) == blocks.Energy(dt), "block curl changed complete-step field evolution");
	}
	const auto finalA = reference.DeviceData(), finalB = blocks.DeviceData();
	Require(finalA.electric == finalB.electric && finalA.magnetic == finalB.magnetic,
		"block fields differ despite matching energy");
	auto copied = blocks;
	auto moved = std::move(copied);
	blocks.Step(dt);
	moved.Step(dt);
	Require(blocks.DeviceData().electric == moved.DeviceData().electric &&
		blocks.DeviceData().magnetic == moved.DeviceData().magnetic, "grid copy/move retained borrowed metric pointers");
	const auto info = blocks.StructuredCurlStorage();
	std::cout << "STRUCTURED_CURL " << name << " rows=" << info.regularRows << " exceptions=" << info.exceptionRows
		<< " bytes=" << info.bytes << std::endl;
}
int main()
{
	try
	{
		IndependentFaceOracle();
		for (const std::string name : {"identity", "local", "graded", "multiple", "material", "uniform-fine", "mixed16", "many"})
		{
			Compare(name);
		}
		std::cout << "PASS: spatial curl blocks, exact sparse exceptions, independent face matrix and copied metrics" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

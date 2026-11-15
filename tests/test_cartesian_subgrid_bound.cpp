/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_inverse_bound.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>

using Graph = std::vector<std::vector<std::size_t>>;

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

static void Check(const Graph &graph)
{
	const std::size_t n = graph.size();
	std::vector<double> factor(n * n, 0);
	for (std::size_t i = 0; i < n; ++i)
	{
		factor[i * n + i] = 1;
		for (std::size_t j : graph[i])
		{
			// Symmetric, nonconstant positive weights; strictly diagonally
			// dominant for these graphs, so the independently built S is SPD.
			factor[i * n + j] = .05 + .03 * ((i + j) % 5);
		}
	}
	for (std::size_t i = 0; i < n; ++i)
	{
		for (std::size_t j = 0; j <= i; ++j)
		{
			double value = factor[i * n + j];
			for (std::size_t k = 0; k < j; ++k)
			{
				value -= factor[i * n + k] * factor[j * n + k];
			}
			Require(i != j || value > 0, "independent test matrix is not SPD");
			factor[i * n + j] = i == j ? std::sqrt(value) : value / factor[j * n + j];
		}
	}
	auto solve = [&](std::vector<double> &rhs)
	{
		for (std::size_t i = 0; i < n; ++i)
		{
			for (std::size_t j = 0; j < i; ++j)
			{
				rhs[i] -= factor[i * n + j] * rhs[j];
			}
			rhs[i] /= factor[i * n + i];
		}
		for (std::size_t i = n; i-- > 0;)
		{
			for (std::size_t j = i + 1; j < n; ++j)
			{
				rhs[i] -= factor[j * n + i] * rhs[j];
			}
			rhs[i] /= factor[i * n + i];
		}
	};
	std::size_t solves = 0;
	double actual = 0;
	const bool supported = CartesianSubgridDetail::BipartiteInverseNorm(
		graph,
		[&](std::vector<double> &rhs)
		{
			++solves;
			solve(rhs);
		},
		actual);
	Require(supported && solves == 1, "bipartite bound did not use exactly one solve");
	std::vector<double> rows(n, 0), column(n, 0);
	for (std::size_t j = 0; j < n; ++j)
	{
		std::fill(column.begin(), column.end(), 0);
		column[j] = 1;
		solve(column);
		for (std::size_t i = 0; i < n; ++i)
		{
			rows[i] += std::abs(column[i]);
		}
	}
	const double expected = *std::max_element(rows.begin(), rows.end());
	Require(std::abs(actual / expected - 1) < 2e-13, "one-RHS norm differs from the full dense inverse");
}

int main()
{
	try
	{
		Check(Graph(17));
		Graph chain(127);
		for (std::size_t i = 1; i < chain.size(); ++i)
		{
			chain[i - 1].push_back(i);
			chain[i].push_back(i - 1);
		}
		Check(chain);
		Graph sheet(48);
		auto link = [&](std::size_t a, std::size_t b)
		{
			sheet[a].push_back(b);
			sheet[b].push_back(a);
		};
		for (std::size_t x = 0; x < 6; ++x)
		{
			for (std::size_t y = 0; y < 8; ++y)
			{
				if (x)
				{
					link(x * 8 + y, (x - 1) * 8 + y);
				}
				if (y)
				{
					link(x * 8 + y, x * 8 + y - 1);
				}
			}
		}
		Check(sheet);
		Check({{1}, {0}, {}, {4, 6}, {3, 5}, {4, 6}, {3, 5}});
		unsigned int solves = 0;
		double bound = 0;
		auto count = [&](std::vector<double> &) { ++solves; };
		Require(!CartesianSubgridDetail::BipartiteInverseNorm(Graph{{1, 2}, {0, 2}, {0, 1}}, count, bound) &&
					solves == 0,
				"odd cycle was not rejected before calling the solve");
		bool rejected = false;
		try
		{
			CartesianSubgridDetail::BipartiteInverseNorm(Graph{{1}, {9}}, count, bound);
		}
		catch (const std::invalid_argument &)
		{
			rejected = true;
		}
		Require(rejected, "invalid graph index was accepted");
		std::cout << "signed inverse bound: dense, disconnected, graded chain/sheet and negative controls passed\n";
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

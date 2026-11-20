/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_loss_rows.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <utility>

struct Rows
{
	std::vector<std::uint32_t> start{0}, edge;
	std::vector<float> weight, low;
};
using Row = std::vector<std::pair<std::uint32_t, double>>;

static void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::runtime_error(message);
	}
}
template <class Action> static void Reject(Action action)
{
	bool rejected = false;
	try
	{
		action();
	}
	catch (const std::exception &)
	{
		rejected = true;
	}
	Require(rejected, "invalid loss-row assembly accepted");
}

static void Compare(std::vector<Row> reference, bool sorted)
{
	Rows actual, expected;
	CartesianSubgridLossRows builder(actual, reference.size());
	for (std::size_t row = 0; row < reference.size(); ++row)
	{
		builder.Count(row, static_cast<std::uint32_t>(reference[row].size()));
	}
	builder.Allocate();
	for (std::size_t row = 0; row < reference.size(); ++row)
	{
		for (const auto &term : reference[row])
		{
			const float high = static_cast<float>(term.second);
			builder.Add(row, term.first, high, static_cast<float>(term.second - double(high)));
		}
		if (sorted)
		{
			std::sort(reference[row].begin(), reference[row].end());
		}
		for (const auto &term : reference[row])
		{
			expected.edge.push_back(term.first);
			expected.weight.push_back(static_cast<float>(term.second));
			expected.low.push_back(static_cast<float>(term.second - double(expected.weight.back())));
		}
		expected.start.push_back(static_cast<std::uint32_t>(expected.edge.size()));
	}
	builder.Complete(sorted);
	Require(actual.start == expected.start && actual.edge == expected.edge, "loss row bounds/order changed");
	Require(actual.weight.size() == expected.weight.size() && actual.low.size() == expected.low.size(),
		"loss row coefficient count changed");
	if (!actual.weight.empty())
	{
		Require(std::memcmp(actual.weight.data(), expected.weight.data(), actual.weight.size() * sizeof(float)) == 0 &&
			std::memcmp(actual.low.data(), expected.low.data(), actual.low.size() * sizeof(float)) == 0,
			"loss row high/residual bits changed");
	}
	Reject([&] { builder.Count(0); });
	Reject([&] { builder.Allocate(); });
	Reject([&] { builder.Complete(); });
}

int main()
{
	try
	{
		Compare({}, false);
		Compare(std::vector<Row>(13), true);
		Compare({{}, {{4, .1}, {0, -0.}, {2, -1. / 3}}, {}, {{1, 1e-20}}}, true);
		Compare({{}, {{7, .1}, {7, -.25}, {1, 0.}}, {}}, false);
		std::vector<Row> generated(4097);
		for (std::size_t row = 0; row < generated.size(); ++row)
		{
			for (std::uint32_t column = 0; column < row % 11; ++column)
			{
				generated[row].emplace_back(20 - column, std::sin(double(row + column)));
			}
		}
		Compare(generated, true);
		Compare(generated, false);
		Rows output;
		CartesianSubgridLossRows builder(output, 2);
		Reject([&] { builder.Count(2); });
		Reject([&] { builder.Add(0, 0, 1, 0); });
		Reject([&] { builder.Complete(); });
		builder.Count(0, 2);
		builder.Allocate();
		Reject([&] { builder.Count(0); });
		Reject([&] { builder.Add(1, 0, 1, 0); });
		builder.Add(0, 7, 1, 0);
		Reject([&] { builder.Complete(); });
		builder.Add(0, 7, 2, 0);
		Reject([&] { builder.Add(0, 7, 3, 0); });
		Reject([&] { builder.Complete(true); });
		CartesianSubgridLossRows overflow(output, 2);
		overflow.Count(0, std::numeric_limits<std::uint32_t>::max());
		Reject([&] { overflow.Count(0); });
		overflow.Count(1);
		Reject([&] { overflow.Allocate(); });
		Require(output.edge.empty() && output.weight.empty() && output.low.empty(),
			"overflow allocated coefficient arrays");
		std::cout << "PASS: loss CSR assembly matches ordered nested rows; invalid/count overflow guards retained\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}

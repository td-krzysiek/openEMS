/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_trace_rows.h"
#include <iostream>
#include <stdexcept>

static void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::runtime_error(message);
	}
}

int main()
{
	try
	{
		CartesianSubgridTraceRows packed;
		std::vector<std::vector<CartesianSubgridTraceRows::Term>> reference;
		Require(packed.empty() && packed.size() == 0, "initial trace rows are not empty");
		for (std::size_t row = 0; row < 4097; ++row)
		{
			std::vector<CartesianSubgridTraceRows::Term> terms;
			for (std::size_t term = 0; term < row % 7; ++term)
			{
				// Repeated indices, signed weights and empty rows must retain order.
				terms.emplace_back((row + term) % 11, (term % 2 ? -1. : 1.) / (row + term + 1));
			}
			packed.push_back(terms);
			reference.push_back(std::move(terms));
		}
		const auto compare = [&](const CartesianSubgridTraceRows &actual)
		{
			Require(actual.size() == reference.size(), "trace row count differs");
			for (std::size_t row = 0; row < reference.size(); ++row)
			{
				const auto view = actual[row];
				Require(view.size() == reference[row].size(), "trace term count differs");
				auto expected = reference[row].begin();
				for (const auto &term : view)
				{
					Require(term == *expected++, "trace term or order changed");
				}
				Require(view.empty() == reference[row].empty(), "empty trace changed");
			}
		};
		compare(packed);
		auto copy = packed;
		packed.clear();
		Require(packed.empty() && packed.size() == 0, "clear retained trace rows");
		compare(copy);
		compare(CartesianSubgridTraceRows(copy));
		packed.push_back({{17, -.25}});
		Require(packed.size() == 1 && packed[0].size() == 1 && packed[0].begin()->first == 17,
			"rebuild after clear failed");
		CartesianSubgridTraceRows singleton;
		for (std::size_t row = 0; row < 32768; ++row)
		{
			singleton.push_back({{row, 1.}});
		}
		Require(singleton.StorageBytes() < 21 * singleton.size(),
			"singleton trace payload lost compact paged storage");
		std::cout << "packed trace rows: exact order, empty/page boundaries, copy and rebuild passed\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_curl.h"

#include <iostream>
#include <numeric>

static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

static void Compare(bool repeated)
{
	std::vector<std::size_t> rows{0}, columns;
	std::vector<double> weights;
	for (std::size_t row = 0; row < 1000; ++row)
	{
		if (row % 19)
		{
			const std::size_t base = row * 3;
			for (std::size_t offset : {0, 1, 63, 64})
			{
				columns.push_back(base + offset);
				weights.push_back((offset % 2 ? -1 : 1) *
								 (repeated ? 1. + .1 * (row % 3) : 1. + row / 1001.));
			}
		}
		rows.push_back(columns.size());
	}
	// A one-ULP variant must remain distinct, including when hashes collide.
	weights[0] = std::nextafter(weights[0], 2.);
	CartesianSubgridCurl curl;
	CartesianSubgridCurl::Builder builder(curl, 1000, 4096);
	std::vector<std::pair<std::size_t, double>> rowEntries;
	for (std::size_t row = 0; row < 1000; ++row)
	{
		rowEntries.clear();
		for (std::size_t i = rows[row]; i < rows[row + 1]; ++i)
		{
			rowEntries.push_back({columns[i], weights[i]});
		}
		builder.Append(rowEntries.begin(), rowEntries.end());
	}
	builder.Finish();
	Require(curl.RowCount() == 1000, "curl lost rows");
	if (repeated)
	{
		Require(curl.PatternCount() == 5, "exact repeated/empty/one-ULP patterns were not preserved");
	}
	std::vector<double> electric(4096), magnetic(1000), expected(1000), transpose(4096), actualTranspose(4096);
	for (std::size_t i = 0; i < electric.size(); ++i)
	{
		electric[i] = std::sin(i * .3791);
	}
	for (std::size_t row = 0; row < magnetic.size(); ++row)
	{
		magnetic[row] = std::cos(row * .219);
		std::size_t next = rows[row];
		curl.ForEach(row, [&](std::size_t column, double value)
		{
			Require(column == columns[next] && value == weights[next], "curl row changed coefficients/term order");
			++next;
		});
		Require(next == rows[row + 1], "curl row lost entries");
		for (std::size_t j = rows[row]; j < rows[row + 1]; ++j)
		{
			expected[row] += weights[j] * electric[columns[j]];
			transpose[columns[j]] += weights[j] * magnetic[row];
		}
	}
	curl.ApplyTranspose(magnetic, actualTranspose);
	curl.Apply(electric, magnetic);
	// This is a same-arithmetic representation check, not cross-engine RF equality.
	Require(expected == magnetic && transpose == actualTranspose, "encoded curl changed same-CPU arithmetic");
	curl.ApplyTranspose(std::vector<double>(1000, 1.), actualTranspose);
	std::vector<double> columnSums(4096, 0);
	for (std::size_t row = 0; row < 1000; ++row)
	{
		for (std::size_t i = rows[row]; i < rows[row + 1]; ++i)
		{
			columnSums[columns[i]] += weights[i];
		}
	}
	Require(columnSums == actualTranspose, "streamed transpose changed coefficient accumulation");
	curl.Apply(electric, magnetic);
	Require(expected == magnetic, "streamed curl changed coefficient accumulation");
	Require(!repeated || curl.AllocationBytes() < (rows.size() + columns.size()) * sizeof(std::size_t) +
			weights.size() * sizeof(double), "repeated curl storage did not shrink");
	CartesianSubgridCurl::Builder emptyBuilder(curl, 32, 2);
	std::vector<std::pair<std::size_t, double>> empty;
	for (std::size_t row = 0; row < 32; ++row)
	{
		emptyBuilder.Append(empty.begin(), empty.end());
	}
	emptyBuilder.Finish();
	Require(curl.PatternCount() == 1, "empty streamed rows lost their common stencil");
	std::vector<double> zeros(32, 19.), field(2, 1.);
	curl.Apply(field, zeros);
	Require(std::all_of(zeros.begin(), zeros.end(), [](double value) { return value == 0; }), "empty-row stencil failed");
	bool rejected = false;
	try
	{
		curl.Apply({}, zeros);
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "invalid field dimensions were accepted");
}

int main()
{
	try
	{
		Compare(true);
		Compare(false);
		for (unsigned int invalid = 0; invalid < 3; ++invalid)
		{
			CartesianSubgridCurl curl;
			bool rejected = false;
			try
			{
				CartesianSubgridCurl::Builder builder(curl, 1, 2);
				std::vector<std::pair<std::size_t, double>> row;
				if (invalid == 0)
				{
					row = {{1, 1.}, {0, 1.}};
				}
				else if (invalid == 1)
				{
					row = {{2, 1.}};
				}
				else
				{
					row = {{1, std::numeric_limits<double>::infinity()}};
				}
				builder.Append(row.begin(), row.end());
			}
			catch (const std::invalid_argument &)
			{
				rejected = true;
			}
			Require(rejected, "malformed curl input was accepted");
		}
		bool incomplete = false;
		try
		{
			CartesianSubgridCurl curl;
			CartesianSubgridCurl::Builder builder(curl, 2, 2);
			std::vector<std::pair<std::size_t, double>> row{{0, 1.}};
			builder.Append(row.begin(), row.end());
			builder.Finish();
		}
		catch (const std::invalid_argument &)
		{
			incomplete = true;
		}
		Require(incomplete, "incomplete streamed curl was accepted");
		std::cout << "translated stencils preserve rows, coefficients, transpose and streamed guards\n";
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
	return 0;
}

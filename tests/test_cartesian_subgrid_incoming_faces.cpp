/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_incoming_faces.h"
#include <array>
#include <iostream>
#include <stdexcept>
#include <utility>

using Rows = CartesianSubgridIncomingFaces;
using Index = Rows::Index;

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
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid incoming-face row accepted");
}

static std::vector<Index> Expected(Index row, bool varied)
{
	std::vector<Index> result;
	if (row % 17 == 0)
	{
		return result;
	}
	const Index base = 32 * row, kind = varied ? row % 7 : 0;
	for (Index term = 0; term < 8; ++term)
	{
		result.push_back(base + term * (kind + 2) + (term * term % 2));
	}
	return result;
}

static void Control(Index count, bool varied)
{
	std::vector<Index> starts{0}, faces;
	for (Index row = 0; row < count; ++row)
	{
		const auto values = Expected(row, varied);
		faces.insert(faces.end(), values.begin(), values.end());
		starts.push_back(static_cast<Index>(faces.size()));
	}
	const std::size_t raw = sizeof(Index) * (starts.size() + faces.size());
	Rows rows(starts, faces, 32 * count + 64);
	// The prepared index owns all addresses after its raw CSR is released.
	std::fill(starts.begin(), starts.end(), 0xdeadbeefu);
	std::fill(faces.begin(), faces.end(), 0xdeadbeefu);
	std::vector<Index>().swap(starts);
	std::vector<Index>().swap(faces);
	Rows copied = rows;
	Rows moved = std::move(copied);
	Require(rows.RowCount() == count && moved.RowCount() == count, "incoming row ownership changed");
	for (Index row = 0; row < count; ++row)
	{
		std::vector<Index> actual;
		moved.ForEach(row, [&](Index face) { actual.push_back(face); });
		Require(actual == Expected(row, varied), "incoming order or integer address changed");
	}
	Reject([&] { moved.ForEach(count, [](Index) {}); });
	Require(rows.PatternCount() <= (varied ? 8u : 2u), "translated rows were not shared");
	if (count >= 4096)
	{
		Require(rows.StorageBytes() < raw / 3, "incoming translated index retained expanded face storage");
	}
	std::cout << "INCOMING_ROWS exact_rows=" << count << " varied=" << varied
		<< " patterns=" << rows.PatternCount() << " raw_bytes=" << raw
		<< " packed_bytes=" << rows.StorageBytes() << std::endl;
}

static void Negative()
{
	Reject([] { Rows rows({}, {}, 1); });
	Reject([] { Rows rows({1}, {}, 1); });
	Reject([] { Rows rows({0, 1}, {}, 1); });
	Reject([] { Rows rows({0, 2, 1}, {0}, 1); });
	Reject([] { Rows rows({0, 2}, {1, 1}, 4); });
	Reject([] { Rows rows({0, 2}, {2, 1}, 4); });
	Reject([] { Rows rows({0, 1}, {4}, 4); });
	Rows empty({0}, {}, 0);
	Require(empty.RowCount() == 0 && empty.PatternCount() == 0, "empty row set changed dimensions");
	Reject([&] { empty.ForEach(0, [](Index) {}); });
}

static void IrregularFallback()
{
	const Index count = 4096;
	std::vector<Index> starts{0}, faces;
	starts.reserve(count + 1);
	faces.reserve(8 * count);
	for (Index row = 0; row < count; ++row)
	{
		for (Index term = 0; term < 8; ++term)
		{
			faces.push_back(term * (row + 1));
		}
		starts.push_back(static_cast<Index>(faces.size()));
	}
	const auto capacity = sizeof(Index) * (starts.capacity() + faces.capacity());
	Rows rows(std::move(starts), std::move(faces), 8 * (count + 1));
	Require(!rows.IsCompressed() && rows.PatternCount() == 0,
		"unique incoming rows did not retain the smaller raw representation");
	Require(rows.StorageBytes() == sizeof(Rows) + capacity,
		"raw fallback retained packed arrays or duplicated CSR");
	for (Index row = 0; row < count; ++row)
	{
		Index term = 0;
		rows.ForEach(row, [&](Index face)
		{
			Require(face == term * (row + 1), "raw fallback changed face order");
			++term;
		});
		Require(term == 8, "raw fallback changed row size");
	}
	std::cout << "INCOMING_FALLBACK exact_rows=" << count << " raw_only=1" << std::endl;
}

int main()
{
	try
	{
		for (Index count : {1u, 17u, 4096u, 131072u})
		{
			Control(count, false);
			Control(count, true);
		}
		Negative();
		IrregularFallback();
		std::cout << "PASS: exact translated incoming-face rows, ownership, bounds and storage" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

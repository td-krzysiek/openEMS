/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_packed_edges.h"
#include <array>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <utility>
#include <vector>

struct Edge
{
	int axis;
	std::array<int, 3> lo;
	int size;
};

static bool Same(const Edge &left, const Edge &right)
{
	return left.axis == right.axis && left.lo == right.lo && left.size == right.size;
}

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
		CartesianSubgridPackedEdges<Edge> packed;
		std::vector<Edge> oracle;
		Require(packed.empty() && packed.begin() == packed.end() && packed.AllocationBytes() == 0,
			"empty packed edges allocate or iterate values");
		const std::array<int, 6> positions{0, 1, 16383, 16384, 262144, 524287};
		for (int axis = 0; axis < 3; ++axis)
		{
			for (int length : {1, 2, 4, 8, 16})
			{
				for (int x : positions)
				{
					for (int y : positions)
					{
						for (int z : positions)
						{
							const Edge edge{axis, {x, y, z}, length};
							packed.push_back(edge);
							oracle.push_back(edge);
							Require(Same(packed.at(oracle.size() - 1), edge), "coordinate bit fields collide");
						}
					}
				}
			}
		}
		packed.push_back(oracle[13]);
		oracle.push_back(oracle[13]);
		Require(packed.size() == oracle.size(), "packed insertion/duplicate order changed");
		std::size_t count = 0;
		for (const auto &edge : packed)
		{
			Require(Same(edge, oracle[count++]), "packed iteration changed coordinates");
		}
		Require(count == oracle.size() && std::size_t(std::distance(packed.begin(), packed.end())) == count,
			"packed iterator length differs");
		auto copy = packed;
		const auto iterator = copy.begin();
		packed.push_back({1, {17, 19, 23}, 8});
		Require(copy.size() == oracle.size() && Same(*iterator, oracle.front()), "copy borrows mutable storage");
		auto moved = std::move(copy);
		Require(moved.size() == oracle.size() && Same(moved.at(13), oracle[13]), "move loses coordinates");
		std::size_t rejected = 0;
		for (const Edge invalid : std::vector<Edge>{{-1, {0, 0, 0}, 1}, {3, {0, 0, 0}, 1},
			{0, {-1, 0, 0}, 1}, {1, {0, -1, 0}, 1}, {2, {0, 0, -1}, 1},
			{0, {524288, 0, 0}, 1}, {1, {0, 524288, 0}, 1}, {2, {0, 0, 524288}, 1},
			{0, {0, 0, 0}, 0}, {0, {0, 0, 0}, 3}, {0, {0, 0, 0}, 32}})
		{
			try
			{
				moved.push_back(invalid);
			}
			catch (const std::invalid_argument &)
			{
				++rejected;
			}
			Require(moved.size() == oracle.size(), "rejected insertion mutated packed edges");
		}
		Require(rejected == 11, "packed edge guard did not reject malformed geometry");
		bool bounds = false;
		try
		{
			moved.at(moved.size());
		}
		catch (const std::out_of_range &)
		{
			bounds = true;
		}
		Require(bounds, "packed edge out-of-range read accepted");
		Require(moved.AllocationBytes() < oracle.capacity() * sizeof(Edge), "packed edge addresses did not reduce storage");
		std::cout << "PACKED_EDGES exact_values=" << oracle.size() << " guards=" << rejected
			<< " coordinate_max=524287 grid_guard_max=262144 copy_move=1 bytes=" << moved.AllocationBytes()
			<< " raw_bytes=" << oracle.capacity() * sizeof(Edge) << '\n';
		std::cout << "PASS: exact packed integer edges, canonical iteration, ownership and checked guards\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

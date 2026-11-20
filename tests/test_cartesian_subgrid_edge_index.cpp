/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_edge_index.h"
#include "FDTD/cartesian_subgrid.h"
#include <iostream>
#include <map>
#include <string>
#include <tuple>

struct Edge
{
	int axis;
	std::array<int, 3> lo;
	int size;
	bool operator<(const Edge &other) const
	{
		return std::tie(axis, lo, size) < std::tie(other.axis, other.lo, other.size);
	}
};

static void Require(bool condition, const char *message)
{
	if (!condition)
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
	catch (const std::exception &)
	{
		return;
	}
	throw std::runtime_error("invalid/absent structured edge was accepted");
}

static void Control(const std::string &name, const std::vector<CartesianSubgrid::Region> &regions,
	const CartesianSubgrid::Index &cells, bool uniformFine = false)
{
	CartesianSubgrid grid(cells, {.001, .0011, .0009}, 8.85418781762e-12, 1.2566370614359173e-6, regions, uniformFine);
	std::vector<CartesianSubgridEdgeIndex::Region> boxes;
	if (uniformFine)
	{
		boxes.push_back({{0, 0, 0}, cells, grid.LatticeRatio()});
	}
	else
	{
		for (const auto &region : regions)
		{
			boxes.push_back({region.box.start, region.box.stop, region.ratio});
		}
	}
	CartesianSubgridEdgeIndex candidate(cells, grid.LatticeRatio(), boxes);
	std::map<Edge, std::size_t> oracle;
	std::vector<Edge> visited;
	for (const auto &leaf : grid.Leaves())
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i < 2; ++i)
			{
				for (unsigned int j = 0; j < 2; ++j)
				{
					Edge edge{int(n), leaf.lo, leaf.size[n]};
					edge.lo[(n + 1) % 3] += int(i) * leaf.size[(n + 1) % 3];
					edge.lo[(n + 2) % 3] += int(j) * leaf.size[(n + 2) % 3];
					const std::size_t id = oracle.size();
					Require(candidate.Insert(edge, id) == oracle.emplace(edge, id).second, "insert/alias differs");
					Require(candidate.At(edge) == oracle.at(edge), "edge ID differs");
					visited.push_back(edge);
				}
			}
		}
	}
	for (auto edge = visited.rbegin(); edge != visited.rend(); ++edge)
	{
		Require(candidate.At(*edge) == oracle.at(*edge), "reverse lookup differs");
	}
	Require(candidate.Size() == oracle.size(), "unique edge count differs");
	CartesianSubgridEdgeIndex sparse(cells, grid.LatticeRatio(), boxes);
	std::size_t count = 0;
	for (const auto &entry : oracle)
	{
		if ((count++ % 3) == 0)
		{
			sparse.Insert(entry.first, entry.second);
		}
		else
		{
			Reject([&] { sparse.At(entry.first); });
		}
	}
	Reject([&] { candidate.At(Edge{3, {0, 0, 0}, 1}); });
	Reject([&] { candidate.At(Edge{0, {-1, 0, 0}, 1}); });
	Reject([&] { candidate.At(Edge{0, {99999, 0, 0}, 1}); });
	Reject([&] { candidate.At(Edge{0, {0, 0, 0}, 0}); });
	Reject([&] { candidate.Insert(visited.front(), CartesianSubgridEdgeIndex::absent); });
	std::cout << name << " leaves=" << grid.CellCount() << " unique_edges=" << candidate.Size()
		<< " structured_bytes=" << candidate.StorageBytes() << " map_accounted_bytes="
		<< oracle.size() * (sizeof(Edge) + sizeof(std::size_t) + 4 * sizeof(void *)) << std::endl;
}

int main()
{
	try
	{
		Control("uniform", {{{{2, 2, 2}, {6, 6, 6}}, {1, 1, 1}}}, {10, 10, 10});
		Control("local222", {{{{2, 2, 2}, {6, 6, 6}}, {2, 2, 2}}}, {8, 8, 8});
		Control("local242", {{{{2, 2, 2}, {6, 6, 6}}, {2, 4, 2}}}, {8, 8, 8});
		Control("independent", {{{{2, 2, 2}, {4, 4, 4}}, {4, 1, 2}},
			{{{5, 5, 5}, {7, 7, 7}}, {1, 2, 4}}}, {10, 10, 10});
		Control("ratio16", {{{{2, 2, 2}, {3, 3, 3}}, {16, 16, 16}}}, {6, 6, 6});
		std::vector<CartesianSubgrid::Region> regions;
		for (unsigned int x = 1; x < 8; x += 2)
		{
			for (unsigned int y = 1; y < 8; y += 2)
			{
				for (unsigned int z = 1; z < 8; z += 2)
				{
					regions.push_back({{{x, y, z}, {x + 1, y + 1, z + 1}}, {1, 2, 1}});
				}
			}
		}
		Control("64regions", regions, {10, 10, 10});
		Control("uniformFine222", {{{{2, 2, 2}, {3, 3, 3}}, {2, 2, 2}}}, {6, 6, 6}, true);
		Control("uniformFine412", {{{{2, 2, 2}, {3, 3, 3}}, {4, 1, 2}}}, {6, 6, 6}, true);
		Reject([] { CartesianSubgridEdgeIndex index({0, 1, 1}, {1, 1, 1}, {}); });
		Reject([] { CartesianSubgridEdgeIndex index({2, 2, 2}, {0, 1, 1}, {}); });
		Reject([] { CartesianSubgridEdgeIndex index({2, 2, 2}, {1, 1, 1}, {{{0, 0, 0}, {1, 1, 1}, {2, 1, 1}}}); });
		Reject([] { CartesianSubgridEdgeIndex index({2, 2, 2}, {1, 1, 1}, {{{1, 1, 1}, {3, 2, 2}, {1, 1, 1}}}); });
		Reject([] { CartesianSubgridEdgeIndex index({2, 2, 2}, {1, 1, 1},
			{{{0, 0, 0}, {1, 1, 1}, {1, 1, 1}}, {{0, 0, 0}, {1, 1, 1}, {1, 1, 1}}}); });
		Reject([] { CartesianSubgridEdgeIndex index({4, 4, 4}, {1, 1, 1},
			{{{0, 0, 0}, {1, 1, 1}, {1, 1, 1}}, {{1, 0, 0}, {2, 1, 1}, {1, 1, 1}}}); });
		std::cout << "PASS: structured edge addressing, aliasing and guards" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

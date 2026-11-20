/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_leaf_index.h"
#include <algorithm>
#include <atomic>
#include <array>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

struct Cell
{
	std::array<int, 3> lo, size;
};
using View = CartesianSubgridLeafIndex<Cell>;
using Index = View::Index;

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
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid leaf address accepted");
}
static bool Same(const Cell &a, const Cell &b)
{
	return a.lo == b.lo && a.size == b.size;
}

static std::size_t Control(const Index &cells, const Index &lattice,
	const std::function<Index(const Index &)> &ratio)
{
	std::vector<Cell> explicitLeaves;
	std::vector<std::size_t> first;
	std::vector<Index> ratios;
	for (unsigned int x = 0; x < cells[0]; ++x)
	{
		for (unsigned int y = 0; y < cells[1]; ++y)
		{
			for (unsigned int z = 0; z < cells[2]; ++z)
			{
				const auto r = ratio({x, y, z});
				first.push_back(explicitLeaves.size());
				ratios.push_back(r);
				const std::array<int, 3> size{int(lattice[0] / r[0]), int(lattice[1] / r[1]), int(lattice[2] / r[2])};
				// Independent retained six-integer reference generator.
				for (unsigned int i = 0; i < r[0]; ++i)
				{
					for (unsigned int j = 0; j < r[1]; ++j)
					{
						for (unsigned int k = 0; k < r[2]; ++k)
						{
							explicitLeaves.push_back({{int(lattice[0] * x) + int(i) * size[0],
								int(lattice[1] * y) + int(j) * size[1], int(lattice[2] * z) + int(k) * size[2]}, size});
						}
					}
				}
			}
		}
	}
	View view(cells, lattice, explicitLeaves.size());
	Reject([&] { view.at(0); });
	Reject([&] { view.First(0); });
	for (unsigned int x = 0; x < cells[0]; ++x)
	{
		for (unsigned int y = 0; y < cells[1]; ++y)
		{
			for (unsigned int z = 0; z < cells[2]; ++z)
			{
				view.AppendRoot({x, y, z}, ratio({x, y, z}));
			}
		}
	}
	view.Finish();
	Require(view.size() == explicitLeaves.size(), "leaf count changed");
	for (std::size_t root = 0; root < first.size(); ++root)
	{
		Require(view.First(root) == first[root] && view.Ratio(root) == ratios[root], "root addressing changed");
	}
	for (std::size_t leaf = 0; leaf < explicitLeaves.size(); ++leaf)
	{
		Require(Same(view[leaf], explicitLeaves[leaf]), "leaf order or integer coordinates changed");
	}
	const std::vector<Cell> iterated(view.begin(), view.end());
	Require(iterated.size() == view.size() && std::equal(iterated.begin(), iterated.end(), explicitLeaves.begin(), Same),
		"range traversal changed canonical order");
	View copied = view;
	View moved = std::move(copied);
	Require(copied.size() == 0 && copied.begin() == copied.end(), "moved-from address range is not empty");
	Reject([&] { copied.at(0); });
	View reassigned;
	reassigned = std::move(moved);
	Require(moved.size() == 0, "move assignment left stale source bounds");
	Reject([&] { moved.at(0); });
	moved = std::move(reassigned);
	view = View();
	for (std::size_t leaf = 0; leaf < explicitLeaves.size(); ++leaf)
	{
		Require(Same(moved.at(leaf), explicitLeaves[leaf]), "copy/move retained borrowed metadata");
	}
	Reject([&] { moved.at(moved.size()); });
	Reject([&] { moved.Ratio(first.size()); });
	Reject([&] { moved.AppendRoot({0, 0, 0}, {1, 1, 1}); });
	if (explicitLeaves.size() == first.size())
	{
		Require(moved.AllocationBytes() == 0, "identity topology unexpectedly allocated per-leaf data");
	}
	return explicitLeaves.size();
}

static void SnapshotCopies()
{
	using Snapshot = CartesianSubgridLeafSnapshot<Cell>;
	for (unsigned int repeat = 0; repeat < 128; ++repeat)
	{
		Snapshot shared;
		std::atomic<bool> start{false}, published{false};
		bool consistent = true;
		std::thread publisher([&]
		{
			while (!start.load())
			{
				std::this_thread::yield();
			}
			Snapshot::Value value = std::make_shared<std::vector<Cell>>(17, Cell{{1, 2, 3}, {4, 8, 16}});
			Snapshot::Value absent;
			published.store(shared.Publish(absent, value));
		});
		start.store(true);
		for (unsigned int i = 0; i < 64; ++i)
		{
			Snapshot copied = shared;
			Snapshot assigned;
			assigned = copied;
			const auto value = assigned.Load();
			if (value)
			{
				consistent = consistent && value->size() == 17 && Same(value->back(), Cell{{1, 2, 3}, {4, 8, 16}});
			}
		}
		publisher.join();
		Require(published.load() && consistent, "concurrent diagnostic copy/publication failed");
		Snapshot copied = shared;
		Snapshot moved(std::move(copied));
		Require(!copied.Load() && moved.Load() == shared.Load(), "moved diagnostic copy lost ownership");
		Snapshot assigned;
		assigned = std::move(moved);
		Require(!moved.Load() && assigned.Load() == shared.Load(), "diagnostic move assignment lost ownership");
		Snapshot::Value absent;
		Require(!assigned.Publish(absent, {}) && absent == shared.Load(),
			"diagnostic copy lost its published immutable owner");
	}
	std::cout << "LEAF_SNAPSHOT concurrent_copy_cases=128" << std::endl;
}

static void Negative()
{
	Reject([] { View bad({0, 4, 4}, {2, 2, 2}, 100); });
	Reject([] { View bad({4, 4, 4}, {3, 2, 2}, 100); });
	Reject([] { View bad({16384, 16384, 16384}, {2, 2, 2}, 100); });
	Reject([] { View bad({4, 4, 4}, {2, 2, 2}, 63); });
	View partial({2, 2, 2}, {16, 8, 4}, 9);
	Reject([&] { partial.AppendRoot({0, 0, 1}, {1, 1, 1}); });
	Reject([&] { partial.AppendRoot({0, 0, 0}, {3, 1, 1}); });
	Reject([&] { partial.AppendRoot({0, 0, 0}, {1, 16, 1}); });
	Reject([&] { partial.AppendRoot({0, 0, 0}, {16, 8, 4}); });
	Require(partial.size() == 0, "rejected append changed topology");
	partial.AppendRoot({0, 0, 0}, {1, 1, 1});
	Reject([&] { partial.AppendRoot({0, 0, 0}, {1, 1, 1}); });
	Reject([&] { partial.Finish(); });
	View identity({1, 1, 1}, {2, 2, 2}, 1);
	Reject([&] { identity.AppendRoot({0, 0, 0}, {2, 1, 1}); });
}

int main()
{
	try
	{
		std::size_t cells = 0;
		unsigned int cases = 0;
		for (unsigned int x : {1u, 2u, 4u, 8u, 16u})
		{
			for (unsigned int y : {1u, 2u, 4u, 8u, 16u})
			{
				for (unsigned int z : {1u, 2u, 4u, 8u, 16u})
				{
					cells += Control({2, 3, 4}, {16, 16, 16}, [=](const Index &) { return Index{x, y, z}; });
					++cases;
				}
			}
		}
		cells += Control({5, 6, 7}, {16, 8, 4}, [](const Index &root)
		{
			return Index{1u << (root[0] % 5), 1u << (root[1] % 4), 1u << (root[2] % 3)};
		});
		++cases;
		Negative();
		SnapshotCopies();
		std::cout << "LEAF_ADDRESS cases=" << cases << " exact_cells=" << cells << std::endl;
		std::cout << "PASS: exact root/child reconstruction, canonical traversal and owning copies" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

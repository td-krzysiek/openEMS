/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_domain_map.h"
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>

using Owners = CartesianSubgridDomainOwners;
using Map = CartesianSubgridDomainMap;
using Index = Owners::Index;

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
	throw std::runtime_error("invalid compact ownership request accepted");
}

static Index Worker(std::size_t state, unsigned int pattern)
{
	return pattern == 0 ? 63 : static_cast<Index>(pattern == 1 ? state / 997 % 8 : state % 64);
}

static void Storage(std::size_t size, unsigned int pattern)
{
	Owners packed(true), dense(false);
	packed.reserve(size);
	dense.reserve(size);
	std::vector<Index> oracle;
	oracle.reserve(size);
	for (std::size_t k = 0; k < size; ++k)
	{
		const auto worker = Worker(k, pattern);
		packed.push_back(worker);
		dense.push_back(worker);
		oracle.push_back(worker);
	}
	Require(packed.size() == size && dense.size() == size, "owner count changed");
	const auto packedBytes = packed.StorageBytes(), denseBytes = dense.StorageBytes();
	if (size >= 10007)
	{
		Require(packedBytes < .6 * denseBytes, "mixed owner tiles lost their bounded storage benefit");
		if (pattern == 0)
		{
			Require(packedBytes < .01 * denseBytes, "uniform ownership allocated per-state labels");
		}
	}
	const auto writes = std::min<std::size_t>(size, 1000);
	for (std::size_t k = 0; k < writes; ++k)
	{
		const auto state = (k * 104729) % size;
		const Index worker = static_cast<Index>((k * 17) % 64);
		packed.Set(state, worker);
		dense.Set(state, worker);
		oracle[state] = worker;
	}
	const auto copy = packed;
	for (std::size_t k = 0; k < size; ++k)
	{
		Require(packed.at(k) == oracle[k] && dense.at(k) == oracle[k] && copy.at(k) == oracle[k],
			"packed/dense/copy owner differs from independent labels");
	}
	if (size)
	{
		packed.Set(size - 1, (oracle.back() + 1) % 64);
		Require(copy.at(size - 1) == oracle.back(), "owner copy shares mutable storage");
	}
	const auto count = packed.size(), bytes = packed.StorageBytes();
	Reject([&] { packed.push_back(64); });
	Reject([&] { packed.Set(0, 64); });
	Reject([&] { packed.Set(size, 0); });
	Reject([&] { packed.at(size); });
	Reject([&] { dense.Set(size, 0); });
	Reject([&] { dense.at(size); });
	Reject([&] { packed.reserve(std::size_t(std::numeric_limits<Index>::max()) + 1); });
	Require(packed.size() == count && packed.StorageBytes() == bytes,
		"rejected ownership request mutated storage");
	std::cout << "OWNER_STORAGE states=" << size << " pattern=" << pattern
		<< " compact_bytes=" << packedBytes << " dense_bytes=" << denseBytes
		<< " exact_labels=1 owning_copy=1 rejected_unchanged=1" << std::endl;
}

static void Addresses()
{
	Owners packed(true);
	std::vector<Index> oracle(10007);
	for (std::size_t k = 0; k < oracle.size(); ++k)
	{
		oracle[k] = static_cast<Index>((17 * k + k / 7) % 3);
		packed.push_back(oracle[k]);
	}
	// Interface and port ownership corrections precede immutable map use.
	for (std::size_t k = 200; k < 800; ++k)
	{
		packed.Set(k, 2);
		oracle[k] = 2;
	}
	std::vector<Index> local(oracle.size(), Map::absent), unused;
	std::vector<std::unique_ptr<Map>> direct, expanded;
	for (Index worker = 0; worker < 3; ++worker)
	{
		direct.emplace_back(std::make_unique<Map>(packed, local, worker, true));
		expanded.emplace_back(std::make_unique<Map>(oracle, unused, worker, false));
		Index next = 0;
		for (bool owned : {true, false})
		{
			for (Index state = 0; state < oracle.size(); ++state)
			{
				if ((oracle[state] == worker) != owned || (!owned && (state + worker) % 11 != 0))
				{
					continue;
				}
				Require(direct.back()->emplace(state, next) && expanded.back()->emplace(state, next),
					"first compact-owner address insertion failed");
				Require(!direct.back()->emplace(state, next + 1), "duplicate compact-owner address replaced a slot");
				++next;
			}
		}
	}
	std::size_t checked = 0;
	for (Index worker = 0; worker < 3; ++worker)
	{
		Require(direct[worker]->size() == expanded[worker]->size(), "compact owner address count differs");
		for (Index state = 0; state < oracle.size(); ++state)
		{
			if (oracle[state] == worker || (state + worker) % 11 == 0)
			{
				Require(direct[worker]->at(state) == expanded[worker]->at(state),
					"compact owner changed local/global numbering");
				++checked;
			}
			else
			{
				Reject([&] { direct[worker]->at(state); });
			}
		}
	}
	std::cout << "OWNER_ADDRESSES states=" << oracle.size() << " exact_lookups=" << checked << std::endl;
}

int main()
{
	try
	{
		for (std::size_t size : {0, 1, 63, 255, 256, 257, 511, 512, 10007, 1000000})
		{
			for (unsigned int pattern = 0; pattern < 3; ++pattern)
			{
				Storage(size, pattern);
			}
		}
		Addresses();
		std::cout << "PASS: tiled worker ownership, dense oracle, corrections, copies and map views" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

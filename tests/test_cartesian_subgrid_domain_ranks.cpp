/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_domain_map.h"
#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_map>

using Owners = CartesianSubgridDomainOwners;
using Ranks = CartesianSubgridDomainRanks;
using Map = CartesianSubgridDomainMap;
using Index = Map::Index;

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
	throw std::runtime_error("invalid owned-rank operation accepted");
}

static Index Worker(std::size_t state, unsigned int pattern)
{
	if (pattern == 0)
	{
		return 63;
	}
	if (pattern == 1)
	{
		return static_cast<Index>(state / 997 % 64);
	}
	if (pattern == 2)
	{
		return static_cast<Index>(state % 64);
	}
	if (pattern == 3)
	{
		return static_cast<Index>((state * 104729 + state / 7) % 64);
	}
	return state == 259 || state == 400 ? 17 : 2;
}

static void RankControl(std::size_t size, unsigned int pattern)
{
	Owners owners(true);
	owners.reserve(size);
	for (std::size_t state = 0; state < size; ++state)
	{
		owners.push_back(Worker(state, pattern));
	}
	const Ranks ranks(owners), copy(ranks);
	// The immutable index owns its finalized labels; a mutable builder no
	// longer participates in map queries or index lifetime.
	if (size)
	{
		owners.Set(size - 1, (Worker(size - 1, pattern) + 1) % 64);
	}
	std::array<std::size_t, 64> counts{};
	for (std::size_t state = 0; state < size; ++state)
	{
		const auto worker = Worker(state, pattern);
		Require(ranks.Owners().at(state) == worker && copy.Owners().at(state) == worker,
			"rank index shares mutable labels");
		Require(ranks.at(state) == counts[worker] && copy.at(state) == counts[worker],
			"rank differs from independent per-owner counters");
		++counts[worker];
	}
	for (Index worker = 0; worker < counts.size(); ++worker)
	{
		Require(ranks.Count(worker) == counts[worker], "owned rank count changed");
	}
	Require(ranks.RankStorageBytes() <= size * sizeof(Index), "rank allocation exceeds dense slots");
	if (size >= 10007)
	{
		Require(ranks.Packed() == (pattern != 2 && pattern != 3), "packed/dense preflight selection changed");
	}
	Reject([&] { ranks.at(size); });
	Reject([&] { ranks.Count(64); });
	std::cout << "OWNED_RANKS states=" << size << " pattern=" << pattern
		<< " packed=" << ranks.Packed() << " rank_bytes=" << ranks.RankStorageBytes()
		<< " dense_bytes=" << size * sizeof(Index) << " exact=1 owning_copy=1" << std::endl;
}

static void MapControl(unsigned int pattern)
{
	Owners owners(true);
	std::vector<Index> labels(10007), slots(labels.size(), Map::absent), unused;
	for (Index state = 0; state < labels.size(); ++state)
	{
		labels[state] = Worker(state, pattern);
		owners.push_back(labels[state]);
	}
	const Ranks ranks(std::move(owners));
	std::size_t lookups = 0, ghosts = 0;
	for (Index worker = 0; worker < 64; ++worker)
	{
		Map ranked(ranks, worker), dense(labels, slots, worker, true), oracle(labels, unused, worker, false);
		std::unordered_map<Index, Index> independent;
		std::vector<Index> owned, borrowed;
		for (Index state = 0; state < labels.size(); ++state)
		{
			if (labels[state] == worker)
			{
				owned.push_back(state);
				Reject([&] { ranked.at(state); });
			}
			else if ((state + worker) % 101 == 0)
			{
				borrowed.push_back(state);
			}
		}
		if (!owned.empty())
		{
			Reject([&] { ranked.emplace(owned.front(), 1); });
			if (owned.size() > 1)
			{
				Reject([&] { ranked.emplace(owned[1], 0); });
			}
			if (!borrowed.empty())
			{
				Reject([&] { ranked.emplace(borrowed.front(), 0); });
			}
			Require(ranked.size() == 0 && ranked.HashedEntries() == 0,
				"rejected insertion changed rank map");
		}
		Index next = 0;
		for (const auto *states : {&owned, &borrowed})
		{
			if (states == &borrowed && !borrowed.empty())
			{
				const auto before = ranked.size();
				Reject([&] { ranked.emplace(borrowed.front(), next + 1); });
				Require(ranked.size() == before && ranked.HashedEntries() == 0,
					"skipped ghost insertion changed rank map");
			}
			for (Index global : *states)
			{
				Require(ranked.emplace(global, next) && dense.emplace(global, next) && oracle.emplace(global, next),
					"first canonical insertion failed");
				independent.emplace(global, next++);
				Require(!ranked.emplace(global, next) && !dense.emplace(global, next) && !oracle.emplace(global, next),
					"duplicate insertion changed an address");
			}
		}
		Require(ranked.size() == independent.size() && ranked.HashedEntries() == borrowed.size(),
			"owned addresses still consume hash entries");
		ghosts += borrowed.size();
		for (Index global = 0; global < labels.size(); ++global)
		{
			const auto found = independent.find(global);
			if (found == independent.end())
			{
				Reject([&] { ranked.at(global); });
			}
			else
			{
				Require(ranked.at(global) == found->second && dense.at(global) == found->second &&
					oracle.at(global) == found->second, "rank/dense/hash local numbering changed");
				++lookups;
			}
			if (labels[global] == worker)
			{
				Require(ranked.Owned().at(global) == independent.at(global), "owned rank view changed");
			}
			else
			{
				Reject([&] { ranked.Owned().at(global); });
			}
		}
		Reject([&] { ranked.emplace(0, Map::absent); });
		Reject([&] { ranked.emplace(labels.size(), 0); });
		Reject([&] { ranked.at(labels.size()); });
	}
	Reject([&] { Map invalid(ranks, 64); });
	std::cout << "OWNED_RANK_MAP pattern=" << pattern << " exact_lookups=" << lookups
		<< " ghost_hash_entries=" << ghosts << " guards_unchanged=1" << std::endl;
}

int main()
{
	try
	{
		for (const std::size_t size : {0, 1, 63, 255, 256, 257, 511, 512, 10007, 1000000})
		{
			for (unsigned int pattern = 0; pattern < 5; ++pattern)
			{
				RankControl(size, pattern);
			}
		}
		for (unsigned int pattern = 0; pattern < 5; ++pattern)
		{
			MapControl(pattern);
		}
		std::cout << "PASS: immutable owned ranks match dense/hash maps and reject noncanonical insertion" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

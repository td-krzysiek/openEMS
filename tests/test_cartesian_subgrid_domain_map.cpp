/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_domain_map.h"
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unordered_map>

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
	throw std::runtime_error("invalid direct-address lookup accepted");
}

static void Control()
{
	std::vector<Index> owners(10007);
	for (Index i = 0; i < owners.size(); ++i)
	{
		owners[i] = (17 * i + i / 7) % 3;
	}
	std::vector<Index> local(owners.size(), Map::absent), unused;
	std::vector<std::unique_ptr<Map>> direct, oracle;
	std::vector<std::unordered_map<Index, Index>> independent(3);
	for (Index worker = 0; worker < 3; ++worker)
	{
		direct.emplace_back(std::make_unique<Map>(owners, local, worker, true));
		oracle.emplace_back(std::make_unique<Map>(owners, unused, worker, false));
		// Querying an owned slot before assignment must not return zero.
		for (Index global = 0; global < owners.size(); ++global)
		{
			if (owners[global] == worker)
			{
				Reject([&] { direct.back()->at(global); });
				break;
			}
		}
	}
	std::size_t ghosts = 0;
	for (Index worker = 0; worker < 3; ++worker)
	{
		Index next = 0;
		for (bool owned : {true, false})
		{
			for (Index global = 0; global < owners.size(); ++global)
			{
				if ((owners[global] == worker) != owned || (!owned && (global + worker) % 11 != 0))
				{
					continue;
				}
				const Index value = next++;
				Require(direct[worker]->emplace(global, value) && oracle[worker]->emplace(global, value),
					"first address insertion failed");
				independent[worker].emplace(global, value);
				Require(!direct[worker]->emplace(global, value + 1) &&
					!oracle[worker]->emplace(global, value + 1), "duplicate insertion replaced an address");
				ghosts += !owned;
			}
		}
	}
	std::size_t directHashed = 0, oracleHashed = 0;
	for (Index worker = 0; worker < 3; ++worker)
	{
		Require(direct[worker]->size() == independent[worker].size() &&
			oracle[worker]->size() == independent[worker].size(), "address map count differs");
		directHashed += direct[worker]->HashedEntries();
		oracleHashed += oracle[worker]->HashedEntries();
		const auto owned = direct[worker]->Owned();
		for (Index global = 0; global < owners.size(); ++global)
		{
			const auto found = independent[worker].find(global);
			if (found == independent[worker].end())
			{
				Reject([&] { direct[worker]->at(global); });
				Reject([&] { oracle[worker]->at(global); });
			}
			else
			{
				Require(direct[worker]->at(global) == found->second &&
					oracle[worker]->at(global) == found->second, "global/local numbering differs");
			}
			if (owners[global] == worker)
			{
				Require(owned.at(global) == independent[worker].at(global), "owned view changed an address");
			}
			else
			{
				Reject([&] { owned.at(global); });
			}
		}
		Reject([&] { direct[worker]->at(static_cast<Index>(owners.size())); });
		Reject([&] { oracle[worker]->at(static_cast<Index>(owners.size())); });
		Reject([&] { direct[worker]->emplace(0, Map::absent); });
		Reject([&] { direct[worker]->emplace(static_cast<Index>(owners.size()), 0); });
	}
	Require(directHashed == ghosts && oracleHashed == owners.size() + ghosts,
		"owned states still require hash entries");
	Reject([&] { Map invalid(owners, unused, 0, true); });
	Reject([&] { Map invalid(owners, local, 64, true); });
	std::cout << "DOMAIN_ADDRESS states=" << owners.size() << " ghosts=" << ghosts
		<< " direct_hashed=" << directHashed << " oracle_hashed=" << oracleHashed << std::endl;
}

int main()
{
	try
	{
		Control();
		std::cout << "PASS: shared owned slots, private ghost indices and strict owned-only views" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

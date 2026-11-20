/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DOMAIN_MAP_H
#define CARTESIAN_SUBGRID_DOMAIN_MAP_H

#include "subgrid_domain_owners.h"
#include "subgrid_domain_ranks.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

//! Setup-only global-to-local lookup; its immutable input index outlives it.
/*! Canonical owned ranks avoid global mutable per-state address slots.
 *  Other workers hash only their own ghost indices. Retained direct-slot and
 *  expanded all-entry hash APIs remain independent controls. None of the
 *  referenced inputs may be moved or replaced while a map is in use.
 */
class CartesianSubgridDomainMap
{
  public:
	using Index = std::uint32_t;
	static constexpr Index absent = std::numeric_limits<Index>::max();

	template <class Owners>
	CartesianSubgridDomainMap(const Owners &owners, std::vector<Index> &local,
		Index worker, bool indexed)
		: m_Owners(owners), m_Local(&local), m_Worker(worker), m_Indexed(indexed)
	{
		if (worker >= 64 || owners.size() > absent || (indexed && local.size() != owners.size()))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid owned-address map dimensions");
		}
	}

	CartesianSubgridDomainMap(const CartesianSubgridDomainRanks &ranks, Index worker)
		: m_Owners(ranks.Owners()), m_Worker(worker), m_Indexed(false), m_Ranks(&ranks)
	{
		if (worker >= 64)
		{
			throw std::invalid_argument("CartesianSubgrid: owned-rank worker is out of range");
		}
	}

	bool emplace(Index global, Index local)
	{
		if (local == absent)
		{
			throw std::invalid_argument("CartesianSubgrid: local address uses the absent sentinel");
		}
		const bool owned = m_Owners.at(global) == m_Worker;
		if (m_Ranks && owned)
		{
			const Index rank = m_Ranks->at(global);
			if (rank < m_AssignedOwned)
			{
				return false;
			}
			if (rank != m_AssignedOwned || local != rank || m_Count != m_AssignedOwned)
			{
				throw std::invalid_argument("CartesianSubgrid: owned rank insertion is not canonical");
			}
			++m_AssignedOwned;
			++m_Count;
			return true;
		}
		if (m_Indexed && owned)
		{
			auto &slot = m_Local->at(global);
			if (slot != absent)
			{
				return false;
			}
			slot = local;
			++m_Count;
			return true;
		}
		if (m_Ranks)
		{
			if (m_Map.find(global) != m_Map.end())
			{
				return false;
			}
			if (m_AssignedOwned != m_Ranks->Count(m_Worker) || local != m_Count)
			{
				throw std::invalid_argument("CartesianSubgrid: ghost rank insertion precedes or skips owned states");
			}
		}
		const bool inserted = m_Map.emplace(global, local).second;
		m_Count += inserted;
		return inserted;
	}

	Index at(Index global) const
	{
		const Index owner = m_Owners.at(global);
		if (m_Ranks && owner == m_Worker)
		{
			const Index rank = m_Ranks->at(global);
			if (rank >= m_AssignedOwned)
			{
				throw std::out_of_range("CartesianSubgrid: unassigned owned rank");
			}
			return rank;
		}
		if (m_Indexed && owner == m_Worker)
		{
			const Index value = m_Local->at(global);
			if (value == absent)
			{
				throw std::out_of_range("CartesianSubgrid: unassigned owned local address");
			}
			return value;
		}
		// The owner lookup validates global range in oracle mode as well.
		return m_Map.at(global);
	}

	Index OwnedAt(Index global) const
	{
		if (m_Owners.at(global) != m_Worker)
		{
			throw std::out_of_range("CartesianSubgrid: constitutive lookup requires an owned state");
		}
		return at(global);
	}

	struct OwnedView
	{
		const CartesianSubgridDomainMap &map;
		Index at(Index global) const
		{
			return map.OwnedAt(global);
		}
	};

	OwnedView Owned() const
	{
		return {*this};
	}

	std::size_t size() const
	{
		return m_Count;
	}

	//! Exact entry count, not an implementation-dependent allocator estimate.
	std::size_t HashedEntries() const
	{
		return m_Map.size();
	}

  private:
	CartesianSubgridDomainOwnerView m_Owners;
	std::vector<Index> *m_Local = nullptr;
	Index m_Worker;
	bool m_Indexed;
	const CartesianSubgridDomainRanks *m_Ranks = nullptr;
	std::size_t m_AssignedOwned = 0;
	std::size_t m_Count = 0;
	std::unordered_map<Index, Index> m_Map;
};

#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DOMAIN_MAP_H
#define CARTESIAN_SUBGRID_DOMAIN_MAP_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

//! Setup-only global-to-local lookup. Owners and shared addresses outlive it.
/*! Each global state has one owning worker, so its owned local index needs
 *  only one shared dense slot. Other workers hash their own ghost indices.
 *  The expanded oracle retains the independent all-entry hash-map path.
 */
class CartesianSubgridDomainMap
{
  public:
	using Index = std::uint32_t;
	static constexpr Index absent = std::numeric_limits<Index>::max();

	CartesianSubgridDomainMap(const std::vector<Index> &owners, std::vector<Index> &local,
		Index worker, bool indexed)
		: m_Owners(owners), m_Local(local), m_Worker(worker), m_Indexed(indexed)
	{
		if (worker >= 64 || owners.size() > absent || (indexed && local.size() != owners.size()))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid owned-address map dimensions");
		}
	}

	bool emplace(Index global, Index local)
	{
		if (local == absent)
		{
			throw std::invalid_argument("CartesianSubgrid: local address uses the absent sentinel");
		}
		const bool owned = m_Owners.at(global) == m_Worker;
		if (m_Indexed && owned)
		{
			auto &slot = m_Local.at(global);
			if (slot != absent)
			{
				return false;
			}
			slot = local;
			++m_Count;
			return true;
		}
		const bool inserted = m_Map.emplace(global, local).second;
		m_Count += inserted;
		return inserted;
	}

	Index at(Index global) const
	{
		const Index owner = m_Owners.at(global);
		if (m_Indexed && owner == m_Worker)
		{
			const Index value = m_Local.at(global);
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
	const std::vector<Index> &m_Owners;
	std::vector<Index> &m_Local;
	Index m_Worker;
	bool m_Indexed;
	std::size_t m_Count = 0;
	std::unordered_map<Index, Index> m_Map;
};

#endif

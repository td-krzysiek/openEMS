/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DOMAIN_RANKS_H
#define CARTESIAN_SUBGRID_DOMAIN_RANKS_H

#include "subgrid_domain_owners.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

//! Immutable ownership and exact ascending-global-ID local owned addresses.
/*! A uniform tile needs only its owner's local base. Mixed tiles keep exact
 *  uint32 ranks. Dense ranks are used when their planned allocation is smaller.
 *  Ownership is moved into the index, so a map cannot pair ranks with labels
 *  from another or subsequently modified ownership assignment. Maps must not
 *  outlive or move their index, just as with the retained direct-address API.
 */
class CartesianSubgridDomainRanks
{
  public:
	using Index = CartesianSubgridDomainOwners::Index;
	static constexpr std::size_t tileSize = CartesianSubgridDomainOwners::tileSize;
	static constexpr std::uint64_t mixedTag = std::uint64_t(1) << 63;

	explicit CartesianSubgridDomainRanks(CartesianSubgridDomainOwners owners)
		: m_Owners(std::move(owners))
	{
		const std::size_t count = m_Owners.size();
		const std::size_t tiles = count / tileSize + (count % tileSize != 0);
		std::size_t mixed = 0;
		for (std::size_t tile = 0; tile < tiles; ++tile)
		{
			mixed += !Uniform(tile);
		}
		const auto denseBytes = Bytes(count, sizeof(Index));
		const auto tileBytes = Bytes(tiles, sizeof(std::uint64_t));
		const auto mixedBytes = Bytes(mixed, sizeof(Tile));
		m_Packed = tileBytes <= denseBytes && mixedBytes < denseBytes - tileBytes;
		if (m_Packed)
		{
			m_Tiles.resize(tiles);
			m_Mixed.reserve(mixed);
			// Account for actual capacities as well as the preflight sizes.
			if (RankStorageBytes() >= denseBytes)
			{
				std::vector<std::uint64_t>().swap(m_Tiles);
				std::vector<Tile>().swap(m_Mixed);
				m_Packed = false;
			}
		}
		if (!m_Packed)
		{
			m_Dense.resize(count);
			for (std::size_t state = 0; state < count; ++state)
			{
				m_Dense[state] = Next(m_Owners.at(state));
			}
			return;
		}
		for (std::size_t tile = 0; tile < tiles; ++tile)
		{
			const auto begin = tile * tileSize;
			const auto end = begin + std::min(tileSize, count - begin);
			if (Uniform(tile))
			{
				const Index worker = m_Owners.at(begin);
				m_Tiles[tile] = m_Counts[worker];
				for (std::size_t state = begin; state < end; ++state)
				{
					Next(worker);
				}
			}
			else
			{
				if (m_Mixed.size() >= mixedTag)
				{
					throw std::invalid_argument("CartesianSubgrid: mixed rank-tile index overflow");
				}
				m_Tiles[tile] = mixedTag | m_Mixed.size();
				Tile values{};
				for (std::size_t state = begin; state < end; ++state)
				{
					values[state - begin] = Next(m_Owners.at(state));
				}
				m_Mixed.push_back(values);
			}
		}
	}

	const CartesianSubgridDomainOwners &Owners() const
	{
		return m_Owners;
	}

	Index at(std::size_t state) const
	{
		if (state >= m_Owners.size())
		{
			throw std::out_of_range("CartesianSubgrid: owned rank state is out of range");
		}
		if (!m_Packed)
		{
			return m_Dense.at(state);
		}
		const auto tile = m_Tiles.at(state / tileSize);
		if (tile & mixedTag)
		{
			return m_Mixed.at(tile & ~mixedTag)[state % tileSize];
		}
		return static_cast<Index>(tile + state % tileSize);
	}

	std::size_t Count(Index worker) const
	{
		return m_Counts.at(worker);
	}

	bool Packed() const
	{
		return m_Packed;
	}

	std::size_t RankStorageBytes() const
	{
		return Bytes(m_Dense.capacity(), sizeof(Index)) + Bytes(m_Tiles.capacity(), sizeof(std::uint64_t)) +
			Bytes(m_Mixed.capacity(), sizeof(Tile));
	}

	std::size_t StorageBytes() const
	{
		return sizeof(*this) + m_Owners.StorageBytes() - sizeof(m_Owners) + RankStorageBytes();
	}

  private:
	using Tile = std::array<Index, tileSize>;
	static std::size_t Bytes(std::size_t count, std::size_t width)
	{
		if (count > std::numeric_limits<std::size_t>::max() / width)
		{
			throw std::invalid_argument("CartesianSubgrid: owned rank storage size overflow");
		}
		return count * width;
	}
	bool Uniform(std::size_t tile) const
	{
		const auto begin = tile * tileSize;
		const auto end = begin + std::min(tileSize, m_Owners.size() - begin);
		const auto owner = m_Owners.at(begin);
		for (auto state = begin + 1; state < end; ++state)
		{
			if (m_Owners.at(state) != owner)
			{
				return false;
			}
		}
		return true;
	}
	Index Next(Index worker)
	{
		auto &count = m_Counts.at(worker);
		if (count >= std::numeric_limits<Index>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: owned rank uses the absent sentinel");
		}
		return static_cast<Index>(count++);
	}
	CartesianSubgridDomainOwners m_Owners;
	std::array<std::size_t, CartesianSubgridDomainOwners::workers> m_Counts{};
	bool m_Packed = false;
	std::vector<Index> m_Dense;
	std::vector<std::uint64_t> m_Tiles;
	std::vector<Tile> m_Mixed;
};

#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DOMAIN_OWNERS_H
#define CARTESIAN_SUBGRID_DOMAIN_OWNERS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

//! Owning setup-only worker labels; field IDs and local addresses stay uint32.
/*! Uniform 256-state tiles store one worker. Mixed tiles store byte labels.
 *  The expanded oracle retains its independent dense uint32 representation.
 *  No geometry approximation, worker reassignment or field renumbering occurs.
 */
class CartesianSubgridDomainOwners
{
  public:
	using Index = std::uint32_t;
	static constexpr Index workers = 64;
	static constexpr std::size_t tileSize = 256;

	explicit CartesianSubgridDomainOwners(bool compact) : m_Compact(compact) {}

	void reserve(std::size_t count)
	{
		if (count > std::numeric_limits<Index>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: owner-label capacity overflow");
		}
		if (m_Compact)
		{
			m_Tiles.reserve(count / tileSize + (count % tileSize != 0));
		}
		else
		{
			m_Dense.reserve(count);
		}
	}

	std::size_t size() const
	{
		return m_Size;
	}

	void push_back(Index worker)
	{
		CheckWorker(worker);
		if (m_Size == std::numeric_limits<Index>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: owner-label count overflow");
		}
		if (!m_Compact)
		{
			m_Dense.push_back(worker);
		}
		else if (m_Size % tileSize == 0)
		{
			m_Tiles.push_back(worker);
		}
		else
		{
			Write(m_Size, worker);
		}
		++m_Size;
	}

	void Set(std::size_t state, Index worker)
	{
		CheckWorker(worker);
		CheckState(state);
		if (m_Compact)
		{
			Write(state, worker);
		}
		else
		{
			m_Dense.at(state) = worker;
		}
	}

	Index at(std::size_t state) const
	{
		CheckState(state);
		if (!m_Compact)
		{
			return m_Dense.at(state);
		}
		const Index tile = m_Tiles.at(state / tileSize);
		return tile < workers ? tile : m_Mixed.at(tile - workers)[state % tileSize];
	}

	Index operator[](std::size_t state) const
	{
		return at(state);
	}

	std::size_t StorageBytes() const
	{
		return sizeof(*this) + m_Dense.capacity() * sizeof(Index) +
			m_Tiles.capacity() * sizeof(Index) + m_Mixed.capacity() * sizeof(Tile);
	}

  private:
	using Tile = std::array<std::uint8_t, tileSize>;
	static void CheckWorker(Index worker)
	{
		if (worker >= workers)
		{
			throw std::invalid_argument("CartesianSubgrid: owner label exceeds64 workers");
		}
	}
	void CheckState(std::size_t state) const
	{
		if (state >= m_Size)
		{
			throw std::out_of_range("CartesianSubgrid: owner state is out of range");
		}
	}
	void Write(std::size_t state, Index worker)
	{
		auto &tile = m_Tiles.at(state / tileSize);
		if (tile < workers)
		{
			if (tile == worker)
			{
				return;
			}
			if (m_Mixed.size() > std::numeric_limits<Index>::max() - workers)
			{
				throw std::invalid_argument("CartesianSubgrid: mixed owner-tile index overflow");
			}
			Tile values;
			values.fill(static_cast<std::uint8_t>(tile));
			m_Mixed.push_back(values);
			tile = workers + static_cast<Index>(m_Mixed.size() - 1);
		}
		m_Mixed.at(tile - workers)[state % tileSize] = static_cast<std::uint8_t>(worker);
	}
	bool m_Compact;
	std::size_t m_Size = 0;
	std::vector<Index> m_Dense, m_Tiles;
	std::vector<Tile> m_Mixed;
};

//! Non-owning reader for retained dense-map tests and owning label storage.
class CartesianSubgridDomainOwnerView
{
  public:
	using Index = std::uint32_t;
	explicit CartesianSubgridDomainOwnerView(const std::vector<Index> &dense) : m_Dense(&dense) {}
	explicit CartesianSubgridDomainOwnerView(const CartesianSubgridDomainOwners &owners) : m_Owners(&owners) {}
	std::size_t size() const
	{
		return m_Owners ? m_Owners->size() : m_Dense->size();
	}
	Index at(std::size_t state) const
	{
		return m_Owners ? m_Owners->at(state) : m_Dense->at(state);
	}

  private:
	const std::vector<Index> *m_Dense = nullptr;
	const CartesianSubgridDomainOwners *m_Owners = nullptr;
};

#endif

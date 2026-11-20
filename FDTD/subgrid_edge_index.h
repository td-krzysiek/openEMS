/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_EDGE_INDEX_H
#define CARTESIAN_SUBGRID_EDGE_INDEX_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

//! Setup/observation address index; no rectangular parent field allocation.
/*! Base-aligned edges share one base address even when seen by a local block.
 *  Other edges belong to a local region's structured address space. A base-cell
 *  owner table resolves interfaces without a hash/tree search. Slots include
 *  boundary/PEC/hanging edges, but only actual independent edges receive IDs.
 */
class CartesianSubgridEdgeIndex
{
  public:
	using Index = std::array<unsigned int, 3>;
	struct Region
	{
		Index start, stop, ratio;
	};
	static constexpr std::uint32_t absent = std::numeric_limits<std::uint32_t>::max();

	CartesianSubgridEdgeIndex(const Index &cells, const Index &ratio, const std::vector<Region> &regions)
		: m_Cells(cells), m_Ratio(ratio)
	{
		std::size_t roots = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (cells[n] == 0 || ratio[n] == 0 || cells[n] > unsigned(std::numeric_limits<int>::max()) / ratio[n])
			{
				throw std::invalid_argument("CartesianSubgrid: invalid structured index extent");
			}
			roots = Product(roots, cells[n]);
		}
		if (regions.size() > 64)
		{
			throw std::invalid_argument("CartesianSubgrid: too many structured index regions");
		}
		m_Owner.assign(roots, -1);
		m_Base = Block({0, 0, 0}, cells, ratio);
		for (const Region &region : regions)
		{
			Index local{}, origin{}, step{};
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (region.start[n] >= region.stop[n] || region.stop[n] > cells[n] ||
					region.ratio[n] == 0 || ratio[n] % region.ratio[n] != 0)
				{
					throw std::invalid_argument("CartesianSubgrid: invalid structured local region");
				}
				local[n] = (region.stop[n] - region.start[n]) * region.ratio[n];
				origin[n] = region.start[n] * ratio[n];
				step[n] = ratio[n] / region.ratio[n];
			}
			for (std::size_t previous = 0; previous < m_Local.size(); ++previous)
			{
				bool separate = false;
				for (unsigned int n = 0; n < 3; ++n)
				{
					separate = separate || region.stop[n] < regions[previous].start[n] ||
						regions[previous].stop[n] < region.start[n];
				}
				if (!separate)
				{
					throw std::invalid_argument("CartesianSubgrid: structured regions need one base-cell separation");
				}
			}
			const auto owner = static_cast<std::int16_t>(m_Local.size());
			m_Local.emplace_back(origin, local, step);
			for (unsigned int x = region.start[0]; x < region.stop[0]; ++x)
			{
				for (unsigned int y = region.start[1]; y < region.stop[1]; ++y)
				{
					for (unsigned int z = region.start[2]; z < region.stop[2]; ++z)
					{
						auto &slot = m_Owner[Root({x, y, z})];
						if (slot != -1)
						{
							throw std::invalid_argument("CartesianSubgrid: overlapping structured index owners");
						}
						slot = owner;
					}
				}
			}
		}
	}

	template <class Edge> bool Insert(const Edge &edge, std::size_t id)
	{
		if (id >= absent)
		{
			throw std::invalid_argument("CartesianSubgrid: structured edge ID exceeds guard");
		}
		const auto address = Locate(edge);
		auto &block = address.block == 0 ? m_Base : m_Local[address.block - 1];
		auto &slot = block.id[edge.axis][address.offset];
		if (slot != absent)
		{
			return false;
		}
		slot = static_cast<std::uint32_t>(id);
		++m_Count;
		return true;
	}

	template <class Edge> std::size_t At(const Edge &edge) const
	{
		const auto address = Locate(edge);
		const auto &block = address.block == 0 ? m_Base : m_Local[address.block - 1];
		const auto value = block.id[edge.axis][address.offset];
		if (value == absent)
		{
			throw std::out_of_range("CartesianSubgrid: absent structured independent edge");
		}
		return value;
	}

	std::size_t Size() const
	{
		return m_Count;
	}

	std::size_t StorageBytes() const
	{
		std::size_t result = m_Owner.capacity() * sizeof(std::int16_t) + m_Base.Bytes() +
			m_Local.capacity() * sizeof(Block);
		for (const auto &block : m_Local)
		{
			result += block.Bytes();
		}
		return result;
	}

  private:
	static std::size_t Product(std::size_t first, std::size_t second)
	{
		if (second != 0 && first > std::numeric_limits<std::size_t>::max() / second)
		{
			throw std::invalid_argument("CartesianSubgrid: structured index allocation overflow");
		}
		return first * second;
	}

	struct Block
	{
		Block() = default;
		Block(const Index &origin, const Index &cells, const Index &step)
			: origin(origin), cells(cells), step(step)
		{
			for (unsigned int n = 0; n < 3; ++n)
			{
				Index shape = cells;
				++shape[(n + 1) % 3];
				++shape[(n + 2) % 3];
				id[n].assign(Product(Product(shape[0], shape[1]), shape[2]), absent);
			}
		}

		template <class Edge> std::optional<std::size_t> Find(const Edge &edge) const
		{
			if (edge.axis < 0 || edge.axis >= 3 || edge.size <= 0 || unsigned(edge.size) != step[edge.axis])
			{
				return std::nullopt;
			}
			Index point{}, shape = cells;
			++shape[(edge.axis + 1) % 3];
			++shape[(edge.axis + 2) % 3];
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (edge.lo[n] < 0 || unsigned(edge.lo[n]) < origin[n])
				{
					return std::nullopt;
				}
				const unsigned int distance = unsigned(edge.lo[n]) - origin[n];
				if (distance % step[n] != 0 || distance / step[n] >= shape[n])
				{
					return std::nullopt;
				}
				point[n] = distance / step[n];
			}
			return (std::size_t(point[0]) * shape[1] + point[1]) * shape[2] + point[2];
		}

		std::size_t Bytes() const
		{
			return (id[0].capacity() + id[1].capacity() + id[2].capacity()) * sizeof(std::uint32_t);
		}

		Index origin{}, cells{}, step{};
		std::array<std::vector<std::uint32_t>, 3> id;
	};

	std::size_t Root(const Index &point) const
	{
		return (std::size_t(point[0]) * m_Cells[1] + point[1]) * m_Cells[2] + point[2];
	}

	struct Address
	{
		std::size_t block, offset;
	};

	template <class Edge> Address Locate(const Edge &edge) const
	{
		if (const auto base = m_Base.Find(edge))
		{
			return {0, *base};
		}
		if (edge.axis < 0 || edge.axis >= 3)
		{
			throw std::out_of_range("CartesianSubgrid: invalid structured edge axis");
		}
		Index cell{};
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (edge.lo[n] < 0)
			{
				throw std::out_of_range("CartesianSubgrid: negative structured edge coordinate");
			}
			cell[n] = unsigned(edge.lo[n]) / m_Ratio[n];
		}
		const unsigned int a = (edge.axis + 1) % 3, b = (edge.axis + 2) % 3;
		for (unsigned int mask = 0; mask < 4; ++mask)
		{
			Index adjacent = cell;
			bool valid = true;
			for (unsigned int k = 0; k < 2; ++k)
			{
				const unsigned int n = k == 0 ? a : b;
				if (mask & (1u << k))
				{
					if (cell[n] == 0 || unsigned(edge.lo[n]) % m_Ratio[n] != 0)
					{
						valid = false;
						break;
					}
					--adjacent[n];
				}
			}
			if (!valid || adjacent[0] >= m_Cells[0] || adjacent[1] >= m_Cells[1] || adjacent[2] >= m_Cells[2])
			{
				continue;
			}
			const auto owner = m_Owner[Root(adjacent)];
			if (owner >= 0)
			{
				if (const auto local = m_Local[owner].Find(edge))
				{
					return {std::size_t(owner) + 1, *local};
				}
			}
		}
		throw std::out_of_range("CartesianSubgrid: edge has no structured address");
	}

	Index m_Cells, m_Ratio;
	Block m_Base;
	std::vector<Block> m_Local;
	std::vector<std::int16_t> m_Owner;
	std::size_t m_Count = 0;
};
#endif

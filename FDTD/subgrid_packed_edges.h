/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_PACKED_EDGES_H
#define CARTESIAN_SUBGRID_PACKED_EDGES_H

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <vector>

//! Exact integer edge geometry; does not compress floating field/material data.
/*! The existing grid guards allow at most 16384 base intervals and refinement
 *  16 per axis, so every lattice coordinate fits in 19 bits. Two axis bits and
 *  three power-of-two length bits fit alongside the three coordinates in 64
 *  bits. Canonical insertion order and duplicates are preserved. Values are
 *  reconstructed by value; no decoded edge reference escapes this container.
 */
template <class Edge> class CartesianSubgridPackedEdges
{
  public:
	void push_back(const Edge &edge)
	{
		if (edge.axis < 0 || edge.axis > 2 || edge.size < 1 || edge.size > 16 ||
			(edge.size & (edge.size - 1)) != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid packed edge axis or length");
		}
		std::uint64_t value = 0;
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (edge.lo[n] < 0 || std::uint64_t(edge.lo[n]) > CoordinateMask)
			{
				throw std::invalid_argument("CartesianSubgrid: packed edge coordinate exceeds guard");
			}
			value |= std::uint64_t(edge.lo[n]) << (19 * n);
		}
		unsigned int exponent = 0;
		for (int length = edge.size; length > 1; length >>= 1)
		{
			++exponent;
		}
		value |= std::uint64_t(edge.axis) << 57;
		value |= std::uint64_t(exponent) << 59;
		m_Values.push_back(value);
	}

	Edge at(std::size_t index) const
	{
		const auto value = m_Values.at(index);
		Edge edge{};
		for (unsigned int n = 0; n < 3; ++n)
		{
			edge.lo[n] = int((value >> (19 * n)) & CoordinateMask);
		}
		edge.axis = int((value >> 57) & 3);
		edge.size = 1 << ((value >> 59) & 7);
		return edge;
	}

	Edge operator[](std::size_t index) const
	{
		return at(index);
	}
	std::size_t size() const
	{
		return m_Values.size();
	}
	bool empty() const
	{
		return m_Values.empty();
	}
	std::size_t AllocationBytes() const
	{
		return m_Values.capacity() * sizeof(std::uint64_t);
	}

	class const_iterator
	{
	  public:
		using iterator_category = std::forward_iterator_tag;
		using value_type = Edge;
		using difference_type = std::ptrdiff_t;
		using pointer = void;
		using reference = Edge;
		const_iterator() = default;
		const_iterator(const CartesianSubgridPackedEdges *owner, std::size_t index)
			: m_Owner(owner), m_Index(index) {}
		Edge operator*() const
		{
			return m_Owner->at(m_Index);
		}
		const_iterator &operator++()
		{
			++m_Index;
			return *this;
		}
		const_iterator operator++(int)
		{
			const auto previous = *this;
			++*this;
			return previous;
		}
		bool operator==(const const_iterator &other) const
		{
			return m_Owner == other.m_Owner && m_Index == other.m_Index;
		}
		bool operator!=(const const_iterator &other) const
		{
			return !(*this == other);
		}
	  private:
		const CartesianSubgridPackedEdges *m_Owner = nullptr;
		std::size_t m_Index = 0;
	};

	const_iterator begin() const
	{
		return {this, 0};
	}
	const_iterator end() const
	{
		return {this, size()};
	}

  private:
	static constexpr std::uint64_t CoordinateMask = (std::uint64_t(1) << 19) - 1;
	std::vector<std::uint64_t> m_Values;
};
#endif

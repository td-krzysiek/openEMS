/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_LEAF_INDEX_H
#define CARTESIAN_SUBGRID_LEAF_INDEX_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

//! Lazy immutable diagnostic cache; copying also uses an atomic snapshot.
/*! A const diagnostic request may publish concurrently with a read-only grid
 *  copy. Plain shared_ptr member copying would race with that publication.
 *  Copies share only immutable coordinates, never the owning address index.
 */
template <class Cell> class CartesianSubgridLeafSnapshot
{
  public:
	using Value = std::shared_ptr<const std::vector<Cell>>;
	CartesianSubgridLeafSnapshot() = default;
	CartesianSubgridLeafSnapshot(const CartesianSubgridLeafSnapshot &other) : m_Value(other.Load()) {}
	CartesianSubgridLeafSnapshot &operator=(const CartesianSubgridLeafSnapshot &other)
	{
		std::atomic_store(&m_Value, other.Load());
		return *this;
	}
	CartesianSubgridLeafSnapshot(CartesianSubgridLeafSnapshot &&other) noexcept
		: m_Value(std::atomic_exchange(&other.m_Value, Value{})) {}
	CartesianSubgridLeafSnapshot &operator=(CartesianSubgridLeafSnapshot &&other) noexcept
	{
		if (this != &other)
		{
			std::atomic_store(&m_Value, std::atomic_exchange(&other.m_Value, Value{}));
		}
		return *this;
	}
	Value Load() const
	{
		return std::atomic_load(&m_Value);
	}
	bool Publish(Value &expected, const Value &value) const
	{
		return std::atomic_compare_exchange_strong(&m_Value, &expected, value);
	}
  private:
	mutable Value m_Value;
};

//! Owning, exact root/child addresses instead of six integers per leaf.
/*! Root and child order are X/Y/Z with Z fastest, matching the reference
 *  generator. Only integer reconstruction changes; floating coordinates,
 *  field ordering and material sampling are untouched. Copies/moves own all
 *  their metadata; no pointers into a surrounding grid are retained.
 */
template <class Cell> class CartesianSubgridLeafIndex
{
  public:
	using Index = std::array<unsigned int, 3>;
	CartesianSubgridLeafIndex() = default;
	CartesianSubgridLeafIndex(const CartesianSubgridLeafIndex &) = default;
	CartesianSubgridLeafIndex &operator=(const CartesianSubgridLeafIndex &) = default;
	CartesianSubgridLeafIndex(CartesianSubgridLeafIndex &&other) noexcept
	{
		Swap(other);
	}
	CartesianSubgridLeafIndex &operator=(CartesianSubgridLeafIndex &&other) noexcept
	{
		if (this != &other)
		{
			CartesianSubgridLeafIndex moved(std::move(other));
			Swap(moved);
		}
		return *this;
	}
	CartesianSubgridLeafIndex(const Index &cells, const Index &lattice, std::size_t expected)
		: m_Cells(cells), m_Lattice(lattice), m_Expected(expected)
	{
		std::size_t roots = 1;
		const auto limit = std::numeric_limits<std::uint32_t>::max();
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (!cells[n] || !Power(lattice[n]) || roots > limit / cells[n] ||
				cells[n] > unsigned(std::numeric_limits<int>::max()) / lattice[n])
			{
				throw std::invalid_argument("CartesianSubgrid: invalid leaf address dimensions");
			}
			roots *= cells[n];
		}
		if (expected < roots || expected > limit)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid planned leaf count");
		}
		m_RootCount = roots;
		m_Identity = expected == roots;
		if (!m_Identity)
		{
			m_Roots.reserve(roots);
			m_LeafRoot.reserve(expected);
		}
	}

	void AppendRoot(const Index &root, const Index &ratio)
	{
		if (m_Finished || m_NextRoot >= m_RootCount || root != DecodeRoot(m_NextRoot))
		{
			throw std::invalid_argument("CartesianSubgrid: leaf roots must follow canonical order");
		}
		unsigned int powers = 0, children = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (!Power(ratio[n]) || ratio[n] > m_Lattice[n])
			{
				throw std::invalid_argument("CartesianSubgrid: invalid local leaf ratio");
			}
			unsigned int exponent = 0;
			for (unsigned int value = ratio[n]; value > 1; value >>= 1)
			{
				++exponent;
			}
			powers |= exponent << (3 * n);
			children *= ratio[n];
		}
		if (children > m_Expected - m_Count || (m_Identity && children != 1))
		{
			throw std::invalid_argument("CartesianSubgrid: generated leaves exceed their plan");
		}
		if (!m_Identity)
		{
			m_Roots.push_back({static_cast<std::uint32_t>(m_Count), static_cast<std::uint16_t>(powers)});
			m_LeafRoot.insert(m_LeafRoot.end(), children, static_cast<std::uint32_t>(m_NextRoot));
		}
		m_Count += children;
		++m_NextRoot;
	}

	void Finish()
	{
		if (m_NextRoot != m_RootCount || m_Count != m_Expected)
		{
			throw std::invalid_argument("CartesianSubgrid: incomplete leaf address plan");
		}
		m_Finished = true;
	}

	Cell at(std::size_t leaf) const
	{
		if (!m_Finished || leaf >= m_Count)
		{
			throw std::out_of_range("CartesianSubgrid: leaf address is unavailable");
		}
		const auto rootID = m_Identity ? leaf : m_LeafRoot[leaf];
		const auto root = DecodeRoot(rootID);
		const unsigned int packed = m_Identity ? 0 : m_Roots[rootID].powers;
		const unsigned int py = (packed >> 3) & 7, pz = (packed >> 6) & 7;
		const auto child = leaf - (m_Identity ? rootID : m_Roots[rootID].first);
		const std::array<std::size_t, 3> offset{child >> (py + pz),
			(child >> pz) & ((std::size_t(1) << py) - 1), child & ((std::size_t(1) << pz) - 1)};
		Cell cell{};
		for (unsigned int n = 0; n < 3; ++n)
		{
			const unsigned int step = m_Lattice[n] >> ((packed >> (3 * n)) & 7);
			cell.lo[n] = static_cast<int>(m_Lattice[n] * root[n] + offset[n] * step);
			cell.size[n] = static_cast<int>(step);
		}
		return cell;
	}

	Cell operator[](std::size_t leaf) const
	{
		return at(leaf);
	}
	std::size_t First(std::size_t root) const
	{
		CheckRoot(root);
		return m_Identity ? root : m_Roots[root].first;
	}
	Index Ratio(std::size_t root) const
	{
		CheckRoot(root);
		const unsigned int packed = m_Identity ? 0 : m_Roots[root].powers;
		return {1u << (packed & 7), 1u << ((packed >> 3) & 7), 1u << ((packed >> 6) & 7)};
	}
	std::size_t size() const
	{
		return m_Count;
	}
	//! Dynamic allocations only; the containing grid accounts for sizeof(*this).
	std::size_t AllocationBytes() const
	{
		return m_Roots.capacity() * sizeof(Root) + m_LeafRoot.capacity() * sizeof(std::uint32_t);
	}

	class ConstIterator
	{
	  public:
		using iterator_category = std::input_iterator_tag;
		using value_type = Cell;
		using difference_type = std::ptrdiff_t;
		using reference = Cell;
		using pointer = void;
		ConstIterator(const CartesianSubgridLeafIndex *owner, std::size_t leaf) : m_Owner(owner), m_Leaf(leaf) {}
		Cell operator*() const
		{
			return m_Owner->at(m_Leaf);
		}
		ConstIterator &operator++()
		{
			++m_Leaf;
			return *this;
		}
		ConstIterator operator++(int)
		{
			const auto previous = *this;
			++*this;
			return previous;
		}
		bool operator==(const ConstIterator &other) const
		{
			return m_Owner == other.m_Owner && m_Leaf == other.m_Leaf;
		}
		bool operator!=(const ConstIterator &other) const
		{
			return !(*this == other);
		}
	  private:
		const CartesianSubgridLeafIndex *m_Owner;
		std::size_t m_Leaf;
	};
	ConstIterator begin() const
	{
		return {this, 0};
	}
	ConstIterator end() const
	{
		return {this, m_Count};
	}

  private:
	struct Root
	{
		std::uint32_t first;
		std::uint16_t powers;
	};
	void Swap(CartesianSubgridLeafIndex &other) noexcept
	{
		using std::swap;
		swap(m_Cells, other.m_Cells);
		swap(m_Lattice, other.m_Lattice);
		swap(m_Expected, other.m_Expected);
		swap(m_RootCount, other.m_RootCount);
		swap(m_NextRoot, other.m_NextRoot);
		swap(m_Count, other.m_Count);
		swap(m_Identity, other.m_Identity);
		swap(m_Finished, other.m_Finished);
		m_Roots.swap(other.m_Roots);
		m_LeafRoot.swap(other.m_LeafRoot);
	}
	static bool Power(unsigned int value)
	{
		return value >= 1 && value <= 16 && (value & (value - 1)) == 0;
	}
	Index DecodeRoot(std::size_t root) const
	{
		return {static_cast<unsigned int>(root / m_Cells[2] / m_Cells[1]),
			static_cast<unsigned int>((root / m_Cells[2]) % m_Cells[1]),
			static_cast<unsigned int>(root % m_Cells[2])};
	}
	void CheckRoot(std::size_t root) const
	{
		if (!m_Finished || root >= m_RootCount)
		{
			throw std::out_of_range("CartesianSubgrid: base leaf address is unavailable");
		}
	}
	Index m_Cells{}, m_Lattice{};
	std::size_t m_Expected = 0, m_RootCount = 0, m_NextRoot = 0, m_Count = 0;
	bool m_Identity = false, m_Finished = false;
	std::vector<Root> m_Roots;
	std::vector<std::uint32_t> m_LeafRoot;
};

#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_INCOMING_FACES_H
#define CARTESIAN_SUBGRID_INCOMING_FACES_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

//! Exact ordered face IDs, sharing only identical translated integer rows.
/*! CSR is released when integer patterns use fewer allocated bytes; otherwise
 *  owned CSR is retained without a duplicate. Dictionaries are setup-only. No weights,
 *  floating-point tolerances or inferred geometric neighbours in this index.
 *  Each row's ascending face order is the canonical transpose sum order.
 */
class CartesianSubgridIncomingFaces
{
  public:
	using Index = std::uint32_t;
	CartesianSubgridIncomingFaces(std::vector<Index> starts,
		std::vector<Index> ids, std::size_t faceCount)
		: m_RawStart(std::move(starts)), m_RawFaces(std::move(ids))
	{
		const auto &start = m_RawStart;
		const auto &faces = m_RawFaces;
		if (start.empty() || start.front() != 0 || start.back() != faces.size() ||
			start.size() - 1 > std::numeric_limits<Index>::max() ||
			faceCount > std::numeric_limits<Index>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: invalid incoming-face dimensions");
		}
		// Validate every raw row even when the allocation gate keeps CSR.
		for (std::size_t row = 0; row + 1 < start.size(); ++row)
		{
			if (start[row + 1] < start[row] || start[row + 1] > faces.size())
			{
				throw std::invalid_argument("CartesianSubgrid: invalid incoming-face row range");
			}
			for (std::size_t k = start[row]; k < start[row + 1]; ++k)
			{
				if (faces[k] >= faceCount || (k > start[row] && faces[k] <= faces[k - 1]))
				{
					throw std::invalid_argument("CartesianSubgrid: incoming faces must be ordered and unique");
				}
			}
		}
		const auto rawIndices = start.capacity() + faces.capacity();
		if (start.size() - 1 >= rawIndices / 2)
		{
			return;
		}
		m_Base.reserve(start.size() - 1);
		m_Pattern.reserve(start.size() - 1);
		m_Start.push_back(0);
		std::unordered_multimap<std::uint64_t, Index> patterns;
		patterns.reserve(std::min<std::size_t>(start.size() - 1, 65536));
		for (std::size_t row = 0; row + 1 < start.size(); ++row)
		{
			const std::size_t first = start[row], stop = start[row + 1];
			const Index base = first == stop ? 0 : faces[first];
			std::uint64_t hash = 14695981039346656037ULL;
			for (std::size_t k = first; k < stop; ++k)
			{
				hash = (hash ^ (faces[k] - base)) * 1099511628211ULL;
			}
			hash ^= stop - first;
			Index pattern = 0;
			bool found = false;
			const auto candidates = patterns.equal_range(hash);
			for (auto match = candidates.first; match != candidates.second; ++match)
			{
				const std::size_t begin = m_Start[match->second];
				if (m_Start[match->second + 1] - begin != stop - first)
				{
					continue;
				}
				std::size_t k = first;
				while (k < stop && faces[k] - base == m_Offset[begin + k - first])
				{
					++k;
				}
				if (k == stop)
				{
					pattern = match->second;
					found = true;
					break;
				}
			}
			if (!found)
			{
				if (stop - first > std::numeric_limits<Index>::max() - m_Offset.size())
				{
					throw std::invalid_argument("CartesianSubgrid: incoming-face pattern overflow");
				}
				pattern = static_cast<Index>(m_Start.size() - 1);
				for (std::size_t k = first; k < stop; ++k)
				{
					m_Offset.push_back(faces[k] - base);
				}
				m_Start.push_back(static_cast<Index>(m_Offset.size()));
				patterns.emplace(hash, pattern);
			}
			m_Base.push_back(base);
			m_Pattern.push_back(pattern);
			if (PackedIndices() >= rawIndices)
			{
				ClearPacked();
				return; // Validated CSR is smaller; keep it without a second copy.
			}
		}
		m_Start.shrink_to_fit();
		m_Offset.shrink_to_fit();
		std::vector<Index>().swap(m_RawStart);
		std::vector<Index>().swap(m_RawFaces);
	}

	template <class Visitor> void ForEach(std::size_t row, Visitor visitor) const
	{
		if (!m_RawStart.empty())
		{
			for (Index k = m_RawStart.at(row); k < m_RawStart.at(row + 1); ++k)
			{
				visitor(m_RawFaces[k]);
			}
			return;
		}
		const auto pattern = m_Pattern.at(row), base = m_Base[row];
		for (Index k = m_Start[pattern]; k < m_Start[pattern + 1]; ++k)
		{
			visitor(base + m_Offset[k]);
		}
	}
	std::size_t RowCount() const
	{
		return m_RawStart.empty() ? m_Base.size() : m_RawStart.size() - 1;
	}
	std::size_t PatternCount() const
	{
		return m_Start.empty() ? 0 : m_Start.size() - 1;
	}
	std::size_t StorageBytes() const
	{
		return sizeof(*this) + sizeof(Index) * (PackedIndices() + m_RawStart.capacity() + m_RawFaces.capacity());
	}

	bool IsCompressed() const
	{
		return m_RawStart.empty();
	}

  private:
	std::size_t PackedIndices() const
	{
		return m_Base.capacity() + m_Pattern.capacity() + m_Start.capacity() + m_Offset.capacity();
	}
	void ClearPacked()
	{
		for (auto *values : {&m_Base, &m_Pattern, &m_Start, &m_Offset})
		{
			std::vector<Index>().swap(*values);
		}
	}
	std::vector<Index> m_Base, m_Pattern, m_Start, m_Offset;
	std::vector<Index> m_RawStart, m_RawFaces;
};

#endif

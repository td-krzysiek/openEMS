/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_CURL_H
#define CARTESIAN_SUBGRID_CURL_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

//! Exact translated stencils for repeated local Cartesian curl rows.
/*!
 * Reuse coefficients and relative edge indices, without rounding coefficients
 * or changing term/row order. Irregular rows may have their own pattern.
 * Rows are streamed without allocating global CSR. Setup-only dictionaries are freed
 * before stepping. This is an intermediate representation, not a complete
 * structured block field layout or an independently qualified GPU backend.
 */
class CartesianSubgridCurl
{
  public:
	//! Stream normalized rows after electric masses are known, without a global CSR.
	class Builder
	{
	  public:
		Builder(CartesianSubgridCurl &target, std::size_t rows, std::size_t edges)
			: m_Target(target), m_Expected(rows)
		{
			if (edges > std::numeric_limits<std::uint32_t>::max() ||
				rows > std::numeric_limits<std::uint32_t>::max())
			{
				throw std::invalid_argument("CartesianSubgrid: streamed curl exceeds 32-bit address guard");
			}
			target.ClearPatterns();
			target.m_Count = 0;
			target.m_EdgeCount = edges;
			target.m_Base.reserve(rows);
			target.m_Pattern.reserve(rows);
			target.m_Start.push_back(0);
			m_Patterns.reserve(std::min<std::size_t>(rows, 65536));
		}

		template <class Iterator> void Append(Iterator first, Iterator last)
		{
			if (m_Target.m_Count >= m_Expected)
			{
				throw std::invalid_argument("CartesianSubgrid: too many streamed curl rows");
			}
			const std::size_t base = first == last ? 0 : first->first;
			std::uint64_t hash = 14695981039346656037ULL;
			std::size_t count = 0, previous = 0;
			for (auto i = first; i != last; ++i)
			{
				if (i->first >= m_Target.m_EdgeCount || !std::isfinite(i->second) ||
					(count && i->first <= previous))
				{
					throw std::invalid_argument("CartesianSubgrid: invalid streamed curl row");
				}
				previous = i->first;
				++count;
				hash = (hash ^ (i->first - base)) * 1099511628211ULL;
				hash = (hash ^ Bits(i->second)) * 1099511628211ULL;
			}
			hash ^= count;
			std::uint32_t pattern = 0;
			bool found = false;
			const auto candidates = m_Patterns.equal_range(hash);
			for (auto match = candidates.first; match != candidates.second; ++match)
			{
				const std::size_t begin = m_Target.m_Start[match->second];
				if (m_Target.m_Start[match->second + 1] - begin != count)
				{
					continue;
				}
				auto entry = first;
				for (std::size_t k = 0; k < count; ++k, ++entry)
				{
					if (entry->first - base != m_Target.m_Offset[begin + k] ||
						Bits(entry->second) != Bits(m_Target.m_Weight[begin + k]))
					{
						break;
					}
				}
				if (entry == last)
				{
					pattern = match->second;
					found = true;
					break;
				}
			}
			if (!found)
			{
				if (count > std::numeric_limits<std::uint32_t>::max() - m_Target.m_Weight.size())
				{
					throw std::invalid_argument("CartesianSubgrid: streamed curl exceeds coefficient guard");
				}
				pattern = static_cast<std::uint32_t>(m_Target.m_Start.size() - 1);
				for (auto i = first; i != last; ++i)
				{
					m_Target.m_Offset.push_back(static_cast<std::uint32_t>(i->first - base));
					m_Target.m_Weight.push_back(i->second);
				}
				m_Target.m_Start.push_back(static_cast<std::uint32_t>(m_Target.m_Weight.size()));
				m_Patterns.emplace(hash, pattern);
			}
			m_Target.m_Base.push_back(static_cast<std::uint32_t>(base));
			m_Target.m_Pattern.push_back(pattern);
			++m_Target.m_Count;
		}

		void Finish()
		{
			if (m_Target.m_Count != m_Expected)
			{
				throw std::invalid_argument("CartesianSubgrid: incomplete streamed curl");
			}
			m_Target.m_Start.shrink_to_fit();
			m_Target.m_Offset.shrink_to_fit();
			m_Target.m_Weight.shrink_to_fit();
			std::unordered_multimap<std::uint64_t, std::uint32_t>().swap(m_Patterns);
		}

	  private:
		CartesianSubgridCurl &m_Target;
		std::size_t m_Expected;
		std::unordered_multimap<std::uint64_t, std::uint32_t> m_Patterns;
	};

	void Apply(const std::vector<double> &electric, std::vector<double> &magnetic) const
	{
		if (electric.size() != m_EdgeCount || magnetic.size() != m_Count)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid curl field dimensions");
		}
		for (std::size_t row = 0; row < m_Count; ++row)
		{
			const std::size_t pattern = m_Pattern[row], base = m_Base[row];
			double value = 0;
			for (std::size_t i = m_Start[pattern]; i < m_Start[pattern + 1]; ++i)
			{
				value += m_Weight[i] * electric[base + m_Offset[i]];
			}
			magnetic[row] = value;
		}
	}

	void ApplyTranspose(const std::vector<double> &magnetic, std::vector<double> &electric) const
	{
		if (electric.size() != m_EdgeCount || magnetic.size() != m_Count)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid transpose field dimensions");
		}
		std::fill(electric.begin(), electric.end(), 0);
		for (std::size_t row = 0; row < m_Count; ++row)
		{
			const std::size_t pattern = m_Pattern[row], base = m_Base[row];
			for (std::size_t i = m_Start[pattern]; i < m_Start[pattern + 1]; ++i)
			{
				electric[base + m_Offset[i]] += m_Weight[i] * magnetic[row];
			}
		}
	}

	template <class Callback> void ForEach(std::size_t row, Callback callback) const
	{
		if (row >= m_Count)
		{
			throw std::out_of_range("CartesianSubgrid: curl row is out of range");
		}
		const std::size_t pattern = m_Pattern[row], base = m_Base[row];
		for (std::size_t i = m_Start[pattern]; i < m_Start[pattern + 1]; ++i)
		{
			callback(base + m_Offset[i], m_Weight[i]);
		}
	}

	std::size_t RowCount() const
	{
		return m_Count;
	}
	std::size_t PatternCount() const
	{
		return m_Start.empty() ? 0 : m_Start.size() - 1;
	}
	std::size_t AllocationBytes() const
	{
		return PatternBytes();
	}

  private:
	static std::uint64_t Bits(double value)
	{
		std::uint64_t bits;
		static_assert(sizeof(bits) == sizeof(value), "64-bit double storage required");
		std::memcpy(&bits, &value, sizeof(bits));
		return bits;
	}
	std::size_t PatternBytes() const
	{
		return (m_Base.capacity() + m_Pattern.capacity() + m_Start.capacity() + m_Offset.capacity()) *
			   sizeof(std::uint32_t) + m_Weight.capacity() * sizeof(double);
	}
	void ClearPatterns()
	{
		std::vector<std::uint32_t>().swap(m_Base);
		std::vector<std::uint32_t>().swap(m_Pattern);
		std::vector<std::uint32_t>().swap(m_Start);
		std::vector<std::uint32_t>().swap(m_Offset);
		std::vector<double>().swap(m_Weight);
	}
	std::size_t m_Count = 0, m_EdgeCount = 0;
	std::vector<std::uint32_t> m_Base, m_Pattern, m_Start, m_Offset;
	std::vector<double> m_Weight;
};
#endif

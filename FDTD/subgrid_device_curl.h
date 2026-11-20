/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_DEVICE_CURL_H
#define SUBGRID_DEVICE_CURL_H

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

//! Exact translated FP32 high/low row patterns for either device curl gather.
/*! No coefficient rounding, reordering or magnetic-history merging. The same
 *  representation applies independently to curl and transpose. Expanded upload
 *  data remains available to the independent arithmetic oracle.
 */
class CartesianSubgridDeviceCurl
{
  public:
	struct Term
	{
		std::uint32_t column;
		float high, low;
	};
	class Builder;
	CartesianSubgridDeviceCurl() = default;
	CartesianSubgridDeviceCurl(const std::vector<std::uint32_t> &rowStart,
		const std::vector<std::uint32_t> &columns, const std::vector<float> &high,
		const std::vector<float> &low, std::size_t inputCount, bool arbitraryColumns = false);

	//! Validate a borrowed package before any backend allocation or dispatch.
	void Validate(std::size_t rows, std::size_t inputs) const
	{
		const auto limit = std::numeric_limits<std::uint32_t>::max();
		if (rows > limit || inputs > limit || base.size() != rows || pattern.size() != rows ||
			start.empty() || start.size() - 1 > limit || offset.size() > limit || start.front() != 0 ||
			start.back() != offset.size() || weight.size() != offset.size() || residual.size() != offset.size() ||
			!std::is_sorted(start.begin(), start.end()))
		{
			throw std::invalid_argument("CartesianSubgrid: malformed compressed domain curl");
		}
		for (std::size_t k = 0; k < weight.size(); ++k)
		{
			if (!NormalOrZero(weight[k]) || !NormalOrZero(residual[k]))
			{
				throw std::invalid_argument("CartesianSubgrid: nonrepresentable compressed curl coefficient");
			}
		}
		for (std::size_t row = 0; row < rows; ++row)
		{
			const std::size_t selected = pattern[row];
			if (selected >= start.size() - 1)
			{
				throw std::invalid_argument("CartesianSubgrid: compressed curl pattern out of range");
			}
			for (std::size_t k = start[selected]; k < start[selected + 1]; ++k)
			{
				if (offset[k] >= inputs || base[row] > inputs - 1 - offset[k])
				{
					throw std::invalid_argument("CartesianSubgrid: compressed curl column out of range");
				}
			}
		}
	}

	std::size_t DeviceBytes() const
	{
		return (start.size() + offset.size() + base.size() + pattern.size()) * sizeof(std::uint32_t) +
			(weight.size() + residual.size()) * sizeof(float);
	}

	std::vector<std::uint32_t> start{0}, offset, base, pattern;
	std::vector<float> weight, residual;

  private:
	static bool NormalOrZero(float value)
	{
		return std::isfinite(value) && (value == 0 || std::abs(value) >= std::numeric_limits<float>::min());
	}

	static std::uint32_t Bits(float value)
	{
		std::uint32_t result;
		static_assert(sizeof(result) == sizeof(value), "32-bit float storage required");
		std::memcpy(&result, &value, sizeof(result));
		return result;
	}
};
//! Setup-only row stream. A completed result has the same exact pattern/order
//! identity as the expanded constructor, without retaining all input terms.
class CartesianSubgridDeviceCurl::Builder
{
  public:
	Builder(std::size_t rows, std::size_t inputCount, bool arbitraryColumns = false)
		: m_Rows(rows), m_InputCount(inputCount), m_Arbitrary(arbitraryColumns)
	{
		const auto limit = std::numeric_limits<std::uint32_t>::max();
		if (rows > limit || inputCount > limit)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid device curl dimensions");
		}
		m_Dictionary.reserve(std::min<std::size_t>(rows + 1, 65536));
		m_Result.base.reserve(rows);
		m_Result.pattern.reserve(rows);
	}

	void Append(const std::vector<Term> &terms)
	{
		if (m_Finished || m_Result.base.size() >= m_Rows ||
			terms.size() > std::numeric_limits<std::uint32_t>::max())
		{
			throw std::invalid_argument("CartesianSubgrid: invalid device curl stream row");
		}
		std::uint32_t origin = terms.empty() ? 0 : terms.front().column;
		if (m_Arbitrary)
		{
			for (const auto &term : terms)
			{
				origin = std::min(origin, term.column);
			}
		}
		std::uint64_t hash = 14695981039346656037ULL;
		for (std::size_t k = 0; k < terms.size(); ++k)
		{
			const auto &term = terms[k];
			if (term.column >= m_InputCount ||
				(!m_Arbitrary && k && term.column <= terms[k - 1].column) ||
				!NormalOrZero(term.high) || !NormalOrZero(term.low))
			{
				throw std::invalid_argument("CartesianSubgrid: invalid device curl term");
			}
			hash = (hash ^ (term.column - origin)) * 1099511628211ULL;
			hash = (hash ^ Bits(term.high)) * 1099511628211ULL;
			hash = (hash ^ Bits(term.low)) * 1099511628211ULL;
		}
		hash ^= terms.size();
		const auto candidates = m_Dictionary.equal_range(hash);
		std::uint32_t selected = 0;
		bool found = false;
		for (auto candidate = candidates.first; candidate != candidates.second; ++candidate)
		{
			const std::size_t begin = m_Result.start[candidate->second];
			if (m_Result.start[candidate->second + 1] - begin != terms.size())
			{
				continue;
			}
			std::size_t k = 0;
			for (; k < terms.size(); ++k)
			{
				const std::size_t entry = begin + k;
				if (m_Result.offset[entry] != terms[k].column - origin ||
					Bits(m_Result.weight[entry]) != Bits(terms[k].high) ||
					Bits(m_Result.residual[entry]) != Bits(terms[k].low))
				{
					break;
				}
			}
			if (k == terms.size())
			{
				selected = candidate->second;
				found = true;
				break;
			}
		}
		if (!found)
		{
			const auto limit = std::numeric_limits<std::uint32_t>::max();
			if (terms.size() > limit - m_Result.offset.size() || m_Result.start.size() - 1 >= limit)
			{
				throw std::invalid_argument("CartesianSubgrid: device curl pattern storage overflows");
			}
			selected = static_cast<std::uint32_t>(m_Result.start.size() - 1);
			for (const auto &term : terms)
			{
				m_Result.offset.push_back(term.column - origin);
				m_Result.weight.push_back(term.high);
				m_Result.residual.push_back(term.low);
			}
			m_Result.start.push_back(static_cast<std::uint32_t>(m_Result.offset.size()));
			m_Dictionary.emplace(hash, selected);
		}
		m_Result.base.push_back(origin);
		m_Result.pattern.push_back(selected);
	}

	CartesianSubgridDeviceCurl Finish()
	{
		if (m_Finished || m_Result.base.size() != m_Rows)
		{
			throw std::invalid_argument("CartesianSubgrid: incomplete or consumed device curl stream");
		}
		m_Finished = true;
		return std::move(m_Result);
	}

  private:
	std::size_t m_Rows, m_InputCount;
	bool m_Arbitrary, m_Finished = false;
	CartesianSubgridDeviceCurl m_Result;
	std::unordered_multimap<std::uint64_t, std::uint32_t> m_Dictionary;
};

inline CartesianSubgridDeviceCurl::CartesianSubgridDeviceCurl(
	const std::vector<std::uint32_t> &rowStart, const std::vector<std::uint32_t> &columns,
	const std::vector<float> &high, const std::vector<float> &low,
	std::size_t inputCount, bool arbitraryColumns)
{
	const auto limit = std::numeric_limits<std::uint32_t>::max();
	if (rowStart.empty() || rowStart.front() != 0 || rowStart.back() != columns.size() ||
		columns.size() > limit || high.size() != columns.size() || low.size() != columns.size() ||
		rowStart.size() - 1 > limit || inputCount > limit || !std::is_sorted(rowStart.begin(), rowStart.end()))
	{
		throw std::invalid_argument("CartesianSubgrid: invalid device curl dimensions");
	}
	Builder stream(rowStart.size() - 1, inputCount, arbitraryColumns);
	std::vector<Term> terms;
	for (std::size_t row = 0; row + 1 < rowStart.size(); ++row)
	{
		terms.clear();
		for (std::size_t k = rowStart[row]; k < rowStart[row + 1]; ++k)
		{
			terms.push_back({columns[k], high[k], low[k]});
		}
		stream.Append(terms);
	}
	*this = stream.Finish();
}

#endif

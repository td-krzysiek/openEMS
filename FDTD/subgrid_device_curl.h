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
#include <vector>

//! Exact translated FP32 high/low row patterns for either device curl gather.
/*! No coefficient rounding, reordering or magnetic-history merging. The same
 *  representation applies independently to curl and transpose. Expanded upload
 *  data remains available to the independent arithmetic oracle.
 */
class CartesianSubgridDeviceCurl
{
  public:
	CartesianSubgridDeviceCurl(const std::vector<std::uint32_t> &rowStart,
		const std::vector<std::uint32_t> &columns, const std::vector<float> &high,
		const std::vector<float> &low, std::size_t inputCount)
	{
		const auto limit = std::numeric_limits<std::uint32_t>::max();
		if (rowStart.empty() || rowStart.front() != 0 || rowStart.back() != columns.size() ||
			columns.size() > limit || high.size() != columns.size() || low.size() != columns.size() ||
			rowStart.size() - 1 > limit || inputCount > limit)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid device curl dimensions");
		}
		std::unordered_multimap<std::uint64_t, std::uint32_t> dictionary;
		dictionary.reserve(std::min<std::size_t>(rowStart.size(), 65536));
		base.reserve(rowStart.size() - 1);
		pattern.reserve(rowStart.size() - 1);
		for (std::size_t row = 0; row + 1 < rowStart.size(); ++row)
		{
			const std::size_t first = rowStart[row], stop = rowStart[row + 1];
			if (first > stop || stop > columns.size())
			{
				throw std::invalid_argument("CartesianSubgrid: invalid device curl row bounds");
			}
			const std::uint32_t origin = first == stop ? 0 : columns[first];
			std::uint64_t hash = 14695981039346656037ULL;
			for (std::size_t k = first; k < stop; ++k)
			{
				if (columns[k] >= inputCount || (k > first && columns[k] <= columns[k - 1]) ||
					!NormalOrZero(high[k]) || !NormalOrZero(low[k]))
				{
					throw std::invalid_argument("CartesianSubgrid: invalid device curl term");
				}
				hash = (hash ^ (columns[k] - origin)) * 1099511628211ULL;
				hash = (hash ^ Bits(high[k])) * 1099511628211ULL;
				hash = (hash ^ Bits(low[k])) * 1099511628211ULL;
			}
			hash ^= stop - first;
			const auto candidates = dictionary.equal_range(hash);
			std::uint32_t selected = 0;
			bool found = false;
			for (auto candidate = candidates.first; candidate != candidates.second; ++candidate)
			{
				const std::size_t begin = start[candidate->second];
				if (start[candidate->second + 1] - begin != stop - first)
				{
					continue;
				}
				std::size_t k = first;
				for (; k < stop; ++k)
				{
					const std::size_t entry = begin + k - first;
					if (offset[entry] != columns[k] - origin || Bits(weight[entry]) != Bits(high[k]) ||
						Bits(residual[entry]) != Bits(low[k]))
					{
						break;
					}
				}
				if (k == stop)
				{
					selected = candidate->second;
					found = true;
					break;
				}
			}
			if (!found)
			{
				selected = static_cast<std::uint32_t>(start.size() - 1);
				for (std::size_t k = first; k < stop; ++k)
				{
					offset.push_back(columns[k] - origin);
					weight.push_back(high[k]);
					residual.push_back(low[k]);
				}
				start.push_back(static_cast<std::uint32_t>(offset.size()));
				dictionary.emplace(hash, selected);
			}
			base.push_back(origin);
			pattern.push_back(selected);
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
#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_ALLOCATION_LIMITS_H
#define SUBGRID_ALLOCATION_LIMITS_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace CartesianSubgridDetail
{
constexpr std::size_t DefaultLeafLimit = 1000000;
//! Twelve local electric edges per leaf bound primary field/state indices.
//! Packed row/pattern offsets retain their separate checked-overflow guards.
//! This representability bound is not a qualified memory or usable-size claim.
constexpr std::size_t MaximumIndexLeafLimit = std::numeric_limits<std::uint32_t>::max() / 12;

inline void ValidateLeafLimit(std::size_t limit)
{
	if (!limit || limit > MaximumIndexLeafLimit)
	{
		throw std::invalid_argument("CartesianSubgrid: leaf limit must fit positive 32-bit local-edge addressing");
	}
}

inline std::size_t LeafProduct(std::size_t count, std::size_t factor, std::size_t limit)
{
	ValidateLeafLimit(limit);
	if (!factor || count > limit / factor)
	{
		throw std::invalid_argument("CartesianSubgrid: configured leaf-cell allocation limit exceeded");
	}
	return count * factor;
}

inline std::size_t RefinedLeafCount(std::size_t leaves, std::size_t selected,
	std::size_t children, std::size_t limit)
{
	ValidateLeafLimit(limit);
	if (leaves > limit || !selected || !children || children - 1 > (limit - leaves) / selected)
	{
		throw std::invalid_argument("CartesianSubgrid: configured leaf-cell allocation limit exceeded");
	}
	return leaves + (children - 1) * selected;
}
} // namespace CartesianSubgridDetail
#endif

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_INTERFACE_LIMITS_H
#define SUBGRID_INTERFACE_LIMITS_H

#include <cstddef>
#include <stdexcept>

namespace CartesianSubgridDetail
{
constexpr std::size_t MaxBaseIntervals = 16384;
constexpr std::size_t SmallInterfaceEdges = 1024;
constexpr std::size_t MaxInterfaceEdges = 16384;
constexpr std::size_t LargeInterfaceHalfBandwidth = 4;

//! Guard before allocating a factor; retain the old bounded general case.
inline std::size_t InterfaceFactorEntries(std::size_t edges, std::size_t bandwidth)
{
	if (!edges || edges > MaxInterfaceEdges || bandwidth >= edges)
	{
		throw std::invalid_argument("CartesianSubgrid: invalid interface size; at most 16384 edges supported");
	}
	if (edges > SmallInterfaceEdges && bandwidth > LargeInterfaceHalfBandwidth)
	{
		throw std::invalid_argument("CartesianSubgrid: large interface requires half-bandwidth at most 4");
	}
	// The validated bounds also make this product representable on 32-bit hosts.
	return edges * (bandwidth + 1);
}

//! Never replace a failed large one-RHS bound with a quadratic column sweep.
inline void RequireInterfaceBoundMethod(std::size_t edges, bool oneSolve)
{
	if (edges > SmallInterfaceEdges && !oneSolve)
	{
		throw std::invalid_argument("CartesianSubgrid: large interface requires a verified one-solve bound");
	}
}
} // namespace CartesianSubgridDetail
#endif

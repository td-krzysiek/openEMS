/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef FDTD_KERNELS_SUBGRID_GEOMETRY_H
#define FDTD_KERNELS_SUBGRID_GEOMETRY_H
#if !defined(__METAL_VERSION__)
#include "subgrid_arithmetic.h"
#endif

//! Spatial addressing is independent of coefficient approximation. The high/low
//! coefficient pairs below are the exact prepared pairs of the sparse oracle.
struct SubgridMetricPair
{
	float high, low;
};
//! Header: blocks, leaves, independent E count, face-pattern offset.
//! Block: cells[3], reserved[3], canonical E-slot offsets[3].
//! Leaf: block, flattened cell slot, six-face mask, exception prefix.
CM_FN bool subgrid_geometry_regular(CM_GLOBAL const uint32_t *geometry, uint32_t face)
{
	const uint32_t leaf = 4 + 9 * geometry[0] + 4 * (face / 6);
	return (geometry[leaf + 2] & (1u << (face % 6))) != 0;
}
struct SubgridFaceLocation
{
	uint32_t block, axis, x, y, z;
};
CM_FN SubgridFaceLocation subgrid_geometry_location(CM_GLOBAL const uint32_t *geometry, uint32_t face)
{
	const uint32_t leaf = 4 + 9 * geometry[0] + 4 * (face / 6);
	const uint32_t block = 4 + 9 * geometry[leaf];
	uint32_t point[3] = {geometry[leaf + 1] / (geometry[block + 1] * geometry[block + 2]),
		(geometry[leaf + 1] / geometry[block + 2]) % geometry[block + 1], geometry[leaf + 1] % geometry[block + 2]};
	const uint32_t axis = (face % 6) / 2;
	point[axis] += face % 2;
	return {block, axis, point[0], point[1], point[2]};
}
CM_FN uint32_t subgrid_geometry_edge_at(CM_GLOBAL const uint32_t *geometry,
	SubgridFaceLocation location, uint32_t term)
{
	uint32_t point[3] = {location.x, location.y, location.z};
	const uint32_t a = (location.axis + 1) % 3, b = (location.axis + 2) % 3;
	const uint32_t direction = term < 2 ? b : a, cross = term < 2 ? a : b;
	point[cross] += (term % 2 == 0 ? 1 : 0);
	const uint32_t sy = geometry[location.block + 1] + (direction == 1 ? 0 : 1);
	const uint32_t sz = geometry[location.block + 2] + (direction == 2 ? 0 : 1);
	return geometry[geometry[location.block + 6 + direction] + (point[0] * sy + point[1]) * sz + point[2]];
}
CM_FN uint32_t subgrid_geometry_edge_index(CM_GLOBAL const uint32_t *geometry,
	uint32_t face, uint32_t term)
{
	return subgrid_geometry_edge_at(geometry, subgrid_geometry_location(geometry, face), term);
}
CM_FN SubgridMetricPair subgrid_geometry_term(CM_GLOBAL const uint32_t *geometry,
	CM_GLOBAL const float *coefficients, uint32_t face, uint32_t term, uint32_t)
{
	const uint32_t pattern = (geometry[geometry[3] + face] >> 8) * 8;
	return {coefficients[pattern + term], coefficients[pattern + 4 + term]};
}
CM_FN uint32_t subgrid_geometry_exception(CM_GLOBAL const uint32_t *geometry, uint32_t face)
{
	const uint32_t leaf = 4 + 9 * geometry[0] + 4 * (face / 6);
	uint32_t row = geometry[leaf + 3];
	for (uint32_t side = 0; side < face % 6; ++side)
	{
		row += (geometry[leaf + 2] & (1u << side)) == 0 ? 1 : 0;
	}
	return row;
}
CM_FN uint32_t subgrid_geometry_edge(CM_GLOBAL const uint32_t *geometry,
	CM_GLOBAL const float *metrics, uint32_t face, uint32_t term,
	CM_PRIVATE SubgridMetricPair &coefficient)
{
	const uint32_t edge = subgrid_geometry_edge_index(geometry, face, term);
	coefficient = edge == 0xffffffffu ? SubgridMetricPair{0, 0} :
		subgrid_geometry_term(geometry, metrics, face, term, edge);
	return edge;
}
CM_FN SubgridMetricPair subgrid_geometry_coefficient(CM_GLOBAL const uint32_t *geometry,
	CM_GLOBAL const float *metrics, uint32_t face, uint32_t edge)
{
	for (uint32_t term = 0; term < 4; ++term)
	{
		if (subgrid_geometry_edge_index(geometry, face, term) == edge)
		{
			return subgrid_geometry_term(geometry, metrics, face, term, edge);
		}
	}
	return {0, 0}; // Host admission proves that every selected transpose term exists.
}
#endif

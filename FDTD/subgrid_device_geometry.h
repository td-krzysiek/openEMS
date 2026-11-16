/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_DEVICE_GEOMETRY_H
#define SUBGRID_DEVICE_GEOMETRY_H

#include "subgrid_regular_curl.h"
#include "gpu_cartesiansubgrid.h"
#include "kernels/subgrid_update.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

//! Setup-only spatial packing with exact prepared coefficient pairs. Regular
//! block rows share coefficient-only patterns; the transpose reads the SAME
//! pair, preserving work pairing without an independent metric approximation.
class CartesianSubgridDeviceGeometry
{
  public:
	std::vector<std::uint32_t> geometry;
	std::vector<float> metrics;
	std::size_t rows = 0;

	CartesianSubgridDeviceGeometry(const CartesianSubgridRegularCurl &regular,
		const CartesianSubgridMetricView &view)
	{
		if (!regular.RowCount())
		{
			return;
		}
		const auto &blocks = regular.Blocks();
		const auto &addresses = regular.Addresses();
		geometry.resize(4 + 9 * blocks.size() + 4 * addresses.size() + view.faces, 0);
		geometry[0] = Narrow(blocks.size());
		geometry[1] = Narrow(addresses.size());
		geometry[2] = Narrow(view.edges);
		geometry[3] = Narrow(4 + 9 * blocks.size() + 4 * addresses.size());
		for (std::size_t b = 0; b < blocks.size(); ++b)
		{
			const auto &block = blocks[b];
			const std::size_t offset = 4 + 9 * b;
			for (unsigned int n = 0; n < 3; ++n)
			{
				geometry[offset + n] = block.cells[n];
				geometry[offset + 6 + n] = Narrow(geometry.size());
				geometry.insert(geometry.end(), block.edge[n].begin(), block.edge[n].end());
			}
		}
		std::map<std::array<std::uint32_t, 8>, std::uint32_t> patterns;
		std::uint32_t exceptions = 0;
		for (std::size_t leaf = 0; leaf < addresses.size(); ++leaf)
		{
			const std::size_t offset = 4 + 9 * blocks.size() + 4 * leaf;
			geometry[offset] = addresses[leaf].block;
			geometry[offset + 1] = addresses[leaf].cell;
			geometry[offset + 2] = regular.Masks()[leaf];
			geometry[offset + 3] = exceptions;
			for (unsigned int face = 0; face < 6; ++face)
			{
				const std::size_t row = 6 * leaf + face;
				if (!regular.IsRegular(row))
				{
					++exceptions;
					continue;
				}
				std::array<float, 8> coefficients{};
				unsigned int count = 0;
				regular.ForEach(row, view, [&](std::size_t edge, double value)
				{
					for (unsigned int term = 0; term < 4; ++term)
					{
						if (subgrid_geometry_edge_index(geometry.data(), Narrow(row), term) != edge)
						{
							continue;
						}
						const float high = static_cast<float>(value), low = static_cast<float>(value - double(high));
						Require(Normal(high) && Normal(low), "nonrepresentable spatial coefficient");
						coefficients[term] = high;
						coefficients[term + 4] = low;
						++count;
					}
				});
				Require(count == 4, "spatial row lost its exact independent edge");
				std::array<std::uint32_t, 8> key;
				std::memcpy(key.data(), coefficients.data(), sizeof(coefficients));
				auto found = patterns.find(key);
				if (found == patterns.end())
				{
					found = patterns.emplace(key, Narrow(metrics.size())).first;
					metrics.insert(metrics.end(), coefficients.begin(), coefficients.end());
				}
				geometry[geometry[3] + row] = found->second;
				++rows;
			}
		}
		Validate(geometry, metrics, view.edges, view.faces);
	}

	static std::size_t Validate(const std::vector<std::uint32_t> &g,
		const std::vector<float> &m, std::size_t edges, std::size_t faces)
	{
		if (g.empty())
		{
			Require(m.empty(), "coefficients without spatial geometry");
			return 0;
		}
		Require(g.size() >= 4 && g[0] > 0 && g[0] <= 65 && g[1] > 0 &&
			g[1] == faces / 6 && faces % 6 == 0 && g[2] == edges,
			"spatial header dimensions disagree");
		const std::size_t head = 4 + 9 * std::size_t(g[0]) + 4 * std::size_t(g[1]) + faces;
		Require(head <= g.size() && g[3] == head - faces && !m.empty() && m.size() % 8 == 0,
			"spatial face-pattern bounds disagree");
		Require(std::all_of(m.begin(), m.end(), Normal), "nonfinite or subnormal spatial coefficient");
		for (std::size_t b = 0; b < g[0]; ++b)
		{
			const std::size_t offset = 4 + 9 * b;
			CartesianSubgridRegularCurl::Index cells{g[offset], g[offset + 1], g[offset + 2]};
			Require(cells[0] && cells[1] && cells[2] && !g[offset + 3] && !g[offset + 4] && !g[offset + 5],
				"invalid spatial block dimensions/reserved fields");
			for (unsigned int n = 0; n < 3; ++n)
			{
				const std::size_t count = CartesianSubgridRegularCurl::CellCount(CartesianSubgridRegularCurl::EdgeShape(cells, n));
				Require(g[offset + 6 + n] >= head && g[offset + 6 + n] <= g.size() && count <= g.size() - g[offset + 6 + n],
					"spatial E slots are out of bounds");
				for (std::size_t i = g[offset + 6 + n]; i < g[offset + 6 + n] + count; ++i)
				{
					Require(g[i] == CartesianSubgridRegularCurl::absent || g[i] < edges, "invalid independent E slot");
				}
			}
		}
		std::size_t rows = 0, exceptions = 0;
		for (std::size_t leaf = 0; leaf < g[1]; ++leaf)
		{
			const std::size_t offset = 4 + 9 * std::size_t(g[0]) + 4 * leaf;
			Require(g[offset] < g[0] && g[offset + 2] < 64 && g[offset + 3] == exceptions,
				"spatial owner/mask/exception prefix is invalid");
			const std::size_t block = 4 + 9 * std::size_t(g[offset]);
			Require(g[offset + 1] < CartesianSubgridRegularCurl::CellCount({g[block], g[block + 1], g[block + 2]}),
				"spatial leaf cell is out of bounds");
			for (unsigned int face = 0; face < 6; ++face)
			{
				if (!(g[offset + 2] & (1u << face)))
				{
					++exceptions;
					continue;
				}
				++rows;
				const std::size_t pattern = g[g[3] + 6 * leaf + face];
				Require(pattern % 8 == 0 && pattern <= m.size() && m.size() - pattern >= 8,
					"spatial coefficient pattern is out of bounds");
				for (unsigned int term = 0; term < 4; ++term)
				{
					Require(subgrid_geometry_edge_index(g.data(), Narrow(6 * leaf + face), term) < edges,
						"selected spatial face has an absent independent E edge");
				}
			}
		}
		return rows;
	}

	//! Check the spatial topology and both exact arithmetic lanes against the
	//! independent expanded upload, including every selected transpose term.
	static void ValidateRows(const CartesianSubgridDeviceData &data)
	{
		if (data.curlGeometry.empty())
		{
			return;
		}
		const auto &g = data.curlGeometry;
		const auto &m = data.curlMetrics;
		for (std::size_t row = 0; row < data.magnetic.size(); ++row)
		{
			if (!subgrid_geometry_regular(g.data(), Narrow(row)))
			{
				continue;
			}
			Require(data.curlStart[row + 1] - data.curlStart[row] == 4, "selected spatial row is not a four-edge curl");
			for (unsigned int term = 0; term < 4; ++term)
			{
				SubgridMetricPair value;
				const auto edge = subgrid_geometry_edge(g.data(), m.data(), Narrow(row), term, value);
				const auto first = data.curlColumn.begin() + data.curlStart[row];
				const auto stop = data.curlColumn.begin() + data.curlStart[row + 1];
				const auto found = std::find(first, stop, edge);
				Require(found != stop, "spatial edge is absent from the expanded curl");
				const std::size_t k = std::size_t(found - data.curlColumn.begin());
				Require(Same(value.high, data.curlWeight[k]) && Same(value.low, data.curlLow[k]),
					"spatial coefficient differs from the exact expanded pair");
			}
		}
		for (std::size_t edge = 0; edge < data.electric.size(); ++edge)
		{
			for (std::size_t k = data.transposeStart[edge]; k < data.transposeStart[edge + 1]; ++k)
			{
				const auto face = data.transposeColumn[k];
				if (!subgrid_geometry_regular(g.data(), face))
				{
					continue;
				}
				unsigned int count = 0;
				for (unsigned int term = 0; term < 4; ++term)
				{
					if (subgrid_geometry_edge_index(g.data(), face, term) == edge)
					{
						++count;
					}
				}
				const auto value = subgrid_geometry_coefficient(g.data(), m.data(), face, Narrow(edge));
				Require(count == 1 && Same(value.high, data.transposeWeight[k]) && Same(value.low, data.transposeLow[k]),
					"spatial transpose is not the exact paired curl");
			}
		}
	}

  private:
	static bool Same(float a, float b)
	{
		return std::memcmp(&a, &b, sizeof(float)) == 0;
	}
	static void Require(bool value, const char *message)
	{
		if (!value)
		{
			throw std::invalid_argument(std::string("CartesianSubgrid: ") + message);
		}
	}
	static std::uint32_t Narrow(std::size_t value)
	{
		Require(value < std::numeric_limits<std::uint32_t>::max(), "spatial device address overflow");
		return static_cast<std::uint32_t>(value);
	}
	static bool Normal(float value)
	{
		return std::isfinite(value) && (value == 0 || std::abs(value) >= std::numeric_limits<float>::min());
	}
};
#endif

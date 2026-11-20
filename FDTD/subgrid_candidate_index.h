/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPENEMS_SUBGRID_CANDIDATE_INDEX_H
#define OPENEMS_SUBGRID_CANDIDATE_INDEX_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

// Conservative preparation-only broad phase. This preserves the input order
// and never replaces the caller's exact geometry/priority tests. Bounds use
// {xmin,xmax,ymin,ymax,zmin,zmax} in the caller's coordinate system.
template <class T> class CartesianSubgridCandidateIndex
{
public:
	using Point = std::array<double, 3>;
	using Bounds = std::array<double, 6>;

	CartesianSubgridCandidateIndex(std::vector<T> objects, const std::vector<Bounds> &bounds,
		const Bounds &domain, double tolerance, unsigned int bins = 16,
		std::size_t maximumEntries = 1048576)
		: m_All(std::move(objects)), m_Domain(domain), m_Bins(bins)
	{
		if (m_All.size() <= 16 || bounds.size() != m_All.size() || bins == 0 || bins > 32 ||
			!std::isfinite(tolerance) || tolerance < 0)
		{
			return;
		}
		double scale = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (!std::isfinite(domain[2 * n]) || !std::isfinite(domain[2 * n + 1]) ||
				!(domain[2 * n + 1] > domain[2 * n]) ||
				!std::isfinite(domain[2 * n + 1] - domain[2 * n]))
			{
				return;
			}
			scale = std::max({scale, std::abs(domain[2 * n]), std::abs(domain[2 * n + 1])});
		}
		for (const Bounds &box : bounds)
		{
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (!std::isfinite(box[2 * n]) || !std::isfinite(box[2 * n + 1]) ||
					box[2 * n] > box[2 * n + 1])
				{
					return;
				}
				scale = std::max({scale, std::abs(box[2 * n]), std::abs(box[2 * n + 1])});
			}
		}
		const double padding = tolerance + 128 * std::numeric_limits<double>::epsilon() * scale;
		if (!std::isfinite(padding))
		{
			return;
		}
		std::vector<std::vector<T>> prepared(std::size_t(bins) * bins * bins);
		std::size_t entries = 0;
		for (std::size_t i = 0; i < bounds.size(); ++i)
		{
			std::array<unsigned int, 3> first{}, last{};
			bool overlaps = true;
			for (unsigned int n = 0; n < 3; ++n)
			{
				const double lo = bounds[i][2 * n] - padding;
				const double hi = bounds[i][2 * n + 1] + padding;
				if (!std::isfinite(lo) || !std::isfinite(hi))
				{
					return;
				}
				overlaps = overlaps && hi >= domain[2 * n] && lo <= domain[2 * n + 1];
				first[n] = Bin(n, std::max(domain[2 * n], std::min(domain[2 * n + 1], lo)));
				last[n] = Bin(n, std::max(domain[2 * n], std::min(domain[2 * n + 1], hi)));
			}
			if (!overlaps)
			{
				continue;
			}
			for (unsigned int x = first[0]; x <= last[0]; ++x)
			{
				for (unsigned int y = first[1]; y <= last[1]; ++y)
				{
					for (unsigned int z = first[2]; z <= last[2]; ++z)
					{
						if (entries == maximumEntries)
						{
							return;
						}
						prepared[Flat({x, y, z})].push_back(m_All[i]);
						++entries;
					}
				}
			}
		}
		m_Buckets = std::move(prepared);
		m_Entries = entries;
	}

	const std::vector<T> &At(const double *point) const
	{
		return Overlapping(point, point);
	}

	const std::vector<T> &Overlapping(const double *lo, const double *hi) const
	{
		if (m_Buckets.empty() || !lo || !hi)
		{
			return m_All;
		}
		std::array<unsigned int, 3> bin{};
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (!std::isfinite(lo[n]) || !std::isfinite(hi[n]) || lo[n] > hi[n] ||
				lo[n] < m_Domain[2 * n] || hi[n] > m_Domain[2 * n + 1])
			{
				return m_All;
			}
			bin[n] = Bin(n, lo[n]);
			if (bin[n] != Bin(n, hi[n]))
			{
				// Avoid query-time merges, allocation and duplicate removal. A
				// multi-bin leaf conservatively reuses the original full list.
				return m_All;
			}
		}
		return m_Buckets[Flat(bin)];
	}

	bool Indexed() const
	{
		return !m_Buckets.empty();
	}

	std::size_t Entries() const
	{
		return m_Entries;
	}

private:
	unsigned int Bin(unsigned int axis, double value) const
	{
		const double fraction = (value - m_Domain[2 * axis]) /
			(m_Domain[2 * axis + 1] - m_Domain[2 * axis]);
		return static_cast<unsigned int>(std::max(0.0,
			std::min(double(m_Bins - 1), std::floor(fraction * m_Bins))));
	}

	std::size_t Flat(const std::array<unsigned int, 3> &bin) const
	{
		return (std::size_t(bin[0]) * m_Bins + bin[1]) * m_Bins + bin[2];
	}

	std::vector<T> m_All;
	Bounds m_Domain;
	unsigned int m_Bins;
	std::vector<std::vector<T>> m_Buckets;
	std::size_t m_Entries = 0;
};

#endif

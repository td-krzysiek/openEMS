/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_LOSS_ROWS_H
#define SUBGRID_LOSS_ROWS_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

//! Setup-only CSR assembly with one cursor per row, not a vector per row.
//! Fill order is retained unless the caller requests unique-column sorting.
template <class Rows> class CartesianSubgridLossRows
{
  public:
	CartesianSubgridLossRows(Rows &output, std::size_t count) : m_Output(output), m_Count(count)
	{
		if (count > std::numeric_limits<std::uint32_t>::max())
		{
			throw std::length_error("CartesianSubgrid loss row count overflows device indexing");
		}
		m_Output.start.assign(count + 1, 0);
		m_Output.edge.clear();
		m_Output.weight.clear();
		m_Output.low.clear();
	}
	void Count(std::size_t row, std::uint32_t count = 1)
	{
		if (m_Phase != Phase::Count || row >= m_Count)
		{
			throw std::invalid_argument("CartesianSubgrid loss row count phase/index mismatch");
		}
		auto &size = m_Output.start[row + 1];
		if (count > std::numeric_limits<std::uint32_t>::max() - size)
		{
			throw std::length_error("CartesianSubgrid loss row length overflows device indexing");
		}
		size += count;
	}
	void Allocate()
	{
		if (m_Phase != Phase::Count)
		{
			throw std::invalid_argument("CartesianSubgrid loss rows are already allocated");
		}
		std::uint64_t total = 0;
		for (std::size_t row = 0; row < m_Count; ++row)
		{
			total += m_Output.start[row + 1];
			if (total > std::numeric_limits<std::uint32_t>::max())
			{
				throw std::length_error("CartesianSubgrid loss terms overflow device indexing");
			}
			m_Output.start[row + 1] = static_cast<std::uint32_t>(total);
		}
		m_Output.edge.resize(static_cast<std::size_t>(total));
		m_Output.weight.resize(static_cast<std::size_t>(total));
		m_Output.low.resize(static_cast<std::size_t>(total));
		m_Cursor.assign(m_Output.start.begin(), m_Output.start.end() - 1);
		m_Phase = Phase::Fill;
	}
	void Add(std::size_t row, std::uint32_t column, float high, float low)
	{
		if (m_Phase != Phase::Fill || row >= m_Count || m_Cursor[row] >= m_Output.start[row + 1])
		{
			throw std::invalid_argument("CartesianSubgrid loss row fill phase/count mismatch");
		}
		const auto at = m_Cursor[row]++;
		m_Output.edge[at] = column;
		m_Output.weight[at] = high;
		m_Output.low[at] = low;
	}
	void Complete(bool sortedUnique = false)
	{
		if (m_Phase != Phase::Fill)
		{
			throw std::invalid_argument("CartesianSubgrid loss rows are not in the fill phase");
		}
		std::vector<Entry> scratch;
		for (std::size_t row = 0; row < m_Count; ++row)
		{
			const auto first = m_Output.start[row], stop = m_Output.start[row + 1];
			if (m_Cursor[row] != stop)
			{
				throw std::invalid_argument("CartesianSubgrid loss row fill is incomplete");
			}
			if (!sortedUnique || stop - first < 2)
			{
				continue;
			}
			scratch.clear();
			for (auto at = first; at < stop; ++at)
			{
				scratch.push_back({m_Output.edge[at], m_Output.weight[at], m_Output.low[at]});
			}
			std::sort(scratch.begin(), scratch.end(), [](const Entry &a, const Entry &b)
			{
				return a.column < b.column;
			});
			for (std::size_t i = 0; i < scratch.size(); ++i)
			{
				if (i && scratch[i - 1].column == scratch[i].column)
				{
					throw std::invalid_argument("CartesianSubgrid conductance columns are not unique");
				}
				m_Output.edge[first + i] = scratch[i].column;
				m_Output.weight[first + i] = scratch[i].high;
				m_Output.low[first + i] = scratch[i].low;
			}
		}
		m_Phase = Phase::Complete;
	}

  private:
	struct Entry
	{
		std::uint32_t column;
		float high, low;
	};
	enum class Phase { Count, Fill, Complete };
	Rows &m_Output;
	std::size_t m_Count;
	std::vector<std::uint32_t> m_Cursor;
	Phase m_Phase = Phase::Count;
};
#endif

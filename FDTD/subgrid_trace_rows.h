/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_TRACE_ROWS_H
#define SUBGRID_TRACE_ROWS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

//! Ordered sparse rows without one heap allocation per physical loss trace.
//! Bounded pages avoid whole-model term/offset copies while appending.
class CartesianSubgridTraceRows
{
  public:
	using Term = std::pair<std::size_t, double>;
	class View
	{
	  public:
		using Iterator = std::vector<Term>::const_iterator;
		View(Iterator first, Iterator last) : m_First(first), m_Last(last) {}
		Iterator begin() const { return m_First; }
		Iterator end() const { return m_Last; }
		std::size_t size() const { return static_cast<std::size_t>(m_Last - m_First); }
		bool empty() const { return m_First == m_Last; }

	  private:
		Iterator m_First, m_Last;
	};

	std::size_t size() const
	{
		return m_Pages.empty() ? 0 : (m_Pages.size() - 1) * RowsPerPage + m_Pages.back().rows;
	}
	bool empty() const { return size() == 0; }
	void clear() { m_Pages.clear(); }
	void push_back(const std::vector<Term> &row)
	{
		if (m_Pages.empty() || m_Pages.back().rows == RowsPerPage)
		{
			m_Pages.emplace_back();
		}
		auto &page = m_Pages.back();
		if (row.size() > std::numeric_limits<std::uint32_t>::max() - page.terms.size())
		{
			throw std::length_error("CartesianSubgrid trace page term count overflows");
		}
		page.terms.insert(page.terms.end(), row.begin(), row.end());
		page.offset[++page.rows] = static_cast<std::uint32_t>(page.terms.size());
	}
	View operator[](std::size_t row) const
	{
		const auto &page = m_Pages[row / RowsPerPage];
		const auto local = row % RowsPerPage;
		return {page.terms.begin() + page.offset[local], page.terms.begin() + page.offset[local + 1]};
	}
	//! Count terms without visiting individual rows.
	std::size_t TermCount() const
	{
		std::size_t count = 0;
		for (const auto &page : m_Pages)
		{
			count += page.terms.size();
		}
		return count;
	}
	//! Accounted page payload; excludes allocator and deque-map bookkeeping.
	std::size_t StorageBytes() const
	{
		std::size_t bytes = m_Pages.size() * sizeof(Page);
		for (const auto &page : m_Pages)
		{
			bytes += page.terms.capacity() * sizeof(Term);
		}
		return bytes;
	}

  private:
	static constexpr std::size_t RowsPerPage = 1024;
	struct Page
	{
		Page() { terms.reserve(RowsPerPage); }
		std::array<std::uint32_t, RowsPerPage + 1> offset{};
		std::vector<Term> terms;
		std::size_t rows = 0;
	};
	std::deque<Page> m_Pages;
};
#endif

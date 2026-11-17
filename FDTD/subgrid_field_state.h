/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_FIELD_STATE_H
#define CARTESIAN_SUBGRID_FIELD_STATE_H

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

//! An exact initial +0 field has its full logical size but no host value array.
//! Point reads never allocate. CPU writes/explicit full imports materialize it.
class CartesianSubgridFieldState
{
  public:
	CartesianSubgridFieldState() = default;
	CartesianSubgridFieldState(const CartesianSubgridFieldState &) = default;
	CartesianSubgridFieldState &operator=(const CartesianSubgridFieldState &) = default;
	CartesianSubgridFieldState(CartesianSubgridFieldState &&other) noexcept
		: m_Count(std::exchange(other.m_Count, 0)), m_Values(std::move(other.m_Values)) {}
	CartesianSubgridFieldState &operator=(CartesianSubgridFieldState &&other) noexcept
	{
		if (this != &other)
		{
			m_Count = std::exchange(other.m_Count, 0);
			m_Values = std::move(other.m_Values);
		}
		return *this;
	}

	void AssignZero(std::size_t count)
	{
		std::vector<double>().swap(m_Values);
		m_Count = count;
	}
	std::size_t size() const { return m_Count; }
	std::size_t capacity() const { return m_Values.capacity(); }
	bool ImplicitZero() const { return m_Values.empty(); }

	double operator[](std::size_t index) const
	{
		return ImplicitZero() ? 0. : m_Values[index];
	}
	double at(std::size_t index) const
	{
		if (index >= m_Count)
		{
			throw std::out_of_range("CartesianSubgrid: field state index is out of range");
		}
		return (*this)[index];
	}
	double &operator[](std::size_t index)
	{
		return MutableDense()[index];
	}
	double &at(std::size_t index)
	{
		if (index >= m_Count)
		{
			throw std::out_of_range("CartesianSubgrid: field state index is out of range");
		}
		return MutableDense()[index];
	}
	std::vector<double> &MutableDense()
	{
		if (m_Values.size() != m_Count)
		{
			m_Values.assign(m_Count, 0.);
		}
		return m_Values;
	}
	const std::vector<double> &Dense() const
	{
		if (m_Values.size() != m_Count)
		{
			throw std::logic_error("CartesianSubgrid: implicit-zero fields have no contiguous host values");
		}
		return m_Values;
	}

  private:
	std::size_t m_Count = 0;
	std::vector<double> m_Values;
};
#endif

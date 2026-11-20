/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_REGULAR_CURL_H
#define CARTESIAN_SUBGRID_REGULAR_CURL_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

//! Borrowed only for the duration of an operation: copying a grid cannot leave
//! its curl storage referring to another grid's metrics.
struct CartesianSubgridMetricView
{
	const double *mass = nullptr, *reluctance = nullptr;
	std::size_t edges = 0, faces = 0;
	//! Opposite faces of one leaf share a metric, never a field/history state.
	bool pairedReluctance = false;
	double Reluctance(std::size_t face) const
	{
		if (!reluctance || face >= faces || (pairedReluctance && faces % 2))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid borrowed face metric view");
		}
		return reluctance[pairedReluctance ? face / 2 : face];
	}
};

//! Spatial block addressing for ordinary four-edge face curls. Interface,
//! PEC-cut and otherwise unresolved rows remain in the algebraic reference.
class CartesianSubgridRegularCurl
{
  public:
	using Index = std::array<unsigned int, 3>;
	static constexpr std::uint32_t absent = std::numeric_limits<std::uint32_t>::max();
	struct Block
	{
		Index cells{};
		std::array<std::vector<double>, 3> length;
		std::array<std::vector<std::uint32_t>, 3> edge;
	};
	struct Address
	{
		std::uint32_t block = 0, cell = 0;
	};

	static std::size_t Product(std::size_t a, std::size_t b)
	{
		if (b && a > std::numeric_limits<std::uint32_t>::max() / b)
		{
			throw std::invalid_argument("CartesianSubgrid: structured curl address overflow");
		}
		return a * b;
	}
	static std::size_t CellCount(const Index &cells)
	{
		return Product(Product(cells[0], cells[1]), cells[2]);
	}
	static Index EdgeShape(const Index &cells, unsigned int axis)
	{
		if (axis > 2 || cells[(axis + 1) % 3] == absent || cells[(axis + 2) % 3] == absent)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid structured curl edge shape");
		}
		Index shape = cells;
		++shape[(axis + 1) % 3];
		++shape[(axis + 2) % 3];
		return shape;
	}
	static std::size_t Offset(const Index &shape, const Index &point)
	{
		if (point[0] >= shape[0] || point[1] >= shape[1] || point[2] >= shape[2])
		{
			throw std::invalid_argument("CartesianSubgrid: structured curl coordinate outside block");
		}
		return (std::size_t(point[0]) * shape[1] + point[1]) * shape[2] + point[2];
	}

	void Configure(std::vector<Block> blocks, std::vector<Address> addresses, std::size_t edges)
	{
		if (blocks.empty() || addresses.empty() || edges >= absent)
		{
			throw std::invalid_argument("CartesianSubgrid: empty or overflowing structured curl layout");
		}
		for (const auto &block : blocks)
		{
			if (!block.cells[0] || !block.cells[1] || !block.cells[2])
			{
				throw std::invalid_argument("CartesianSubgrid: structured curl block must have volume");
			}
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (block.length[n].size() != block.cells[n] ||
					block.edge[n].size() != CellCount(EdgeShape(block.cells, n)))
				{
					throw std::invalid_argument("CartesianSubgrid: structured curl block dimensions disagree");
				}
				for (double length : block.length[n])
				{
					if (!(length > 0) || !std::isfinite(length))
					{
						throw std::invalid_argument("CartesianSubgrid: structured curl needs positive finite lengths");
					}
				}
				for (std::uint32_t edge : block.edge[n])
				{
					if (edge != absent && edge >= edges)
					{
						throw std::invalid_argument("CartesianSubgrid: structured curl independent edge is out of range");
					}
				}
			}
		}
		for (const Address &address : addresses)
		{
			if (address.block >= blocks.size() || address.cell >= CellCount(blocks[address.block].cells))
			{
				throw std::invalid_argument("CartesianSubgrid: structured curl leaf address is out of range");
			}
		}
		m_Blocks = std::move(blocks);
		m_Address = std::move(addresses);
		m_Mask.assign(m_Address.size(), 0);
		m_Prefix.assign(m_Address.size(), 0);
		m_Count = 0;
	}
	void ResetRows(std::size_t rows)
	{
		if (!m_Blocks.empty() && Product(m_Address.size(), 6) != rows)
		{
			throw std::invalid_argument("CartesianSubgrid: structured curl leaf/face dimensions disagree");
		}
		std::fill(m_Mask.begin(), m_Mask.end(), 0);
		m_Count = 0;
	}
	void BeginRow(std::size_t row, std::size_t exceptions)
	{
		if (!m_Blocks.empty() && row % 6 == 0)
		{
			if (exceptions >= absent)
			{
				throw std::invalid_argument("CartesianSubgrid: too many structured curl exceptions");
			}
			m_Prefix[row / 6] = static_cast<std::uint32_t>(exceptions);
		}
	}
	template <class Iterator>
	bool Admit(std::size_t row, Iterator first, Iterator stop, const CartesianSubgridMetricView &metrics)
	{
		if (m_Blocks.empty())
		{
			return false;
		}
		std::array<std::pair<std::size_t, double>, 4> terms;
		if (!Terms(row, metrics, terms))
		{
			return false;
		}
		for (const auto &term : terms)
		{
			if (first == stop || first->first != term.first ||
				std::memcmp(&first->second, &term.second, sizeof(double)) != 0)
			{
				return false;
			}
			++first;
		}
		if (first != stop)
		{
			return false;
		}
		m_Mask[row / 6] |= static_cast<unsigned char>(1u << (row % 6));
		++m_Count;
		return true;
	}
	bool IsRegular(std::size_t row) const
	{
		return m_Count && (m_Mask.at(row / 6) & (1u << (row % 6)));
	}
	std::size_t ExceptionRow(std::size_t row) const
	{
		if (!m_Count)
		{
			return row;
		}
		std::size_t result = m_Prefix.at(row / 6);
		for (unsigned int side = 0; side < row % 6; ++side)
		{
			result += !(m_Mask[row / 6] & (1u << side));
		}
		return result;
	}
	template <class Callback>
	void ForEach(std::size_t row, const CartesianSubgridMetricView &metrics, Callback callback) const
	{
		std::array<std::pair<std::size_t, double>, 4> terms;
		if (!IsRegular(row) || !Terms(row, metrics, terms))
		{
			throw std::logic_error("CartesianSubgrid: prepared structured curl lost its metrics or topology");
		}
		for (const auto &term : terms)
		{
			callback(term.first, term.second);
		}
	}
	std::size_t RowCount() const
	{
		return m_Count;
	}
	std::size_t BlockCount() const
	{
		return m_Blocks.size();
	}
	std::size_t Bytes() const
	{
		std::size_t bytes = m_Blocks.capacity() * sizeof(Block) + m_Address.capacity() * sizeof(Address) +
			m_Mask.capacity() * sizeof(unsigned char) + m_Prefix.capacity() * sizeof(std::uint32_t);
		for (const auto &block : m_Blocks)
		{
			for (unsigned int n = 0; n < 3; ++n)
			{
				bytes += block.length[n].capacity() * sizeof(double) + block.edge[n].capacity() * sizeof(std::uint32_t);
			}
		}
		return bytes;
	}
	void Finish()
	{
		if (!m_Count)
		{
			std::vector<Block>().swap(m_Blocks);
			std::vector<Address>().swap(m_Address);
			std::vector<unsigned char>().swap(m_Mask);
			std::vector<std::uint32_t>().swap(m_Prefix);
		}
	}

	const std::vector<Block> &Blocks() const
	{
		return m_Blocks;
	}
	const std::vector<Address> &Addresses() const
	{
		return m_Address;
	}
	const std::vector<unsigned char> &Masks() const
	{
		return m_Mask;
	}

  private:
	bool Terms(std::size_t row, const CartesianSubgridMetricView &metrics,
		std::array<std::pair<std::size_t, double>, 4> &terms) const
	{
		if (!metrics.mass || !metrics.reluctance || row >= metrics.faces || row / 6 >= m_Address.size())
		{
			throw std::invalid_argument("CartesianSubgrid: structured curl requires current mass/face metrics");
		}
		const double reluctance = metrics.Reluctance(row);
		if (!(reluctance > 0) || !std::isfinite(reluctance))
		{
			throw std::invalid_argument("CartesianSubgrid: structured curl face metric is invalid");
		}
		const auto &address = m_Address[row / 6];
		const auto &block = m_Blocks[address.block];
		Index cell{address.cell / (block.cells[1] * block.cells[2]),
			(address.cell / block.cells[2]) % block.cells[1], address.cell % block.cells[2]};
		const unsigned int n = unsigned(row % 6) / 2, side = unsigned(row % 2);
		const unsigned int a = (n + 1) % 3, b = (n + 2) % 3;
		const unsigned int direction[4] = {b, b, a, a};
		const unsigned int cross[4] = {a, a, b, b};
		const unsigned int high[4] = {1, 0, 1, 0};
		const double sign[4] = {1, -1, -1, 1};
		for (unsigned int k = 0; k < 4; ++k)
		{
			Index point = cell;
			point[n] += side;
			point[cross[k]] += high[k];
			const std::uint32_t edge = block.edge[direction[k]][Offset(EdgeShape(block.cells, direction[k]), point)];
			if (edge == absent)
			{
				return false;
			}
			if (edge >= metrics.edges || !(metrics.mass[edge] > 0) || !std::isfinite(metrics.mass[edge]))
			{
				throw std::invalid_argument("CartesianSubgrid: structured curl mass is invalid");
			}
			terms[k] = {edge, sign[k] * block.length[direction[k]][cell[direction[k]]] *
				(reluctance / metrics.mass[edge])};
			if (!std::isfinite(terms[k].second))
			{
				throw std::invalid_argument("CartesianSubgrid: structured curl coefficient overflows");
			}
		}
		std::sort(terms.begin(), terms.end(), [](const auto &a, const auto &b)
		{
			return a.first < b.first;
		});
		for (unsigned int k = 1; k < terms.size(); ++k)
		{
			if (terms[k].first <= terms[k - 1].first)
			{
				return false;
			}
		}
		return true;
	}
	std::vector<Block> m_Blocks;
	std::vector<Address> m_Address;
	std::vector<unsigned char> m_Mask;
	std::vector<std::uint32_t> m_Prefix;
	std::size_t m_Count = 0;
};
#endif

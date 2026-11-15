/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/gpu_cartesiansubgrid.h"
#include "FDTD/cartesian_subgrid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace
{
//! Independent FP32 oracle used only by device regression tests.
class Reference_CartesianSubgrid : public GPU_CartesianSubgrid
{
  public:
	explicit Reference_CartesianSubgrid(const CartesianSubgrid &grid)
		: m_Data(grid.DeviceData()), m_Increment(m_Data.electric.size(), 0)
	{
	}

	void Step(double timestep, double sourceAmplitude) override
	{
		const float dt = static_cast<float>(timestep), amplitude = static_cast<float>(sourceAmplitude);
		if (!(timestep > 0) || timestep > m_Data.stableTimestep || !std::isfinite(dt) ||
			dt < std::numeric_limits<float>::min() || !std::isfinite(amplitude))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid reference device step or drive");
		}
		for (std::size_t i = 0; i < m_Increment.size(); ++i)
		{
			float value = 0;
			for (std::size_t k = m_Data.transposeStart[i]; k < m_Data.transposeStart[i + 1]; ++k)
			{
				value += m_Data.transposeWeight[k] * m_Data.magnetic[m_Data.transposeColumn[k]];
			}
			m_Increment[i] = value;
		}
		for (std::size_t block = 0; block < m_Data.bandwidth.size(); ++block)
		{
			const std::size_t first = m_Data.blockStart[block];
			const std::size_t n = m_Data.blockStart[block + 1] - first;
			const std::size_t band = m_Data.bandwidth[block], stride = band + 1, offset = m_Data.lowerStart[block];
			for (std::size_t i = 0; i < n; ++i)
			{
				float value = m_Increment[m_Data.blockEdge[first + i]];
				for (std::size_t j = i > band ? i - band : 0; j < i; ++j)
				{
					value -= m_Data.lower[offset + i * stride + i - j] * m_Increment[m_Data.blockEdge[first + j]];
				}
				m_Increment[m_Data.blockEdge[first + i]] = value / m_Data.lower[offset + i * stride];
			}
			for (std::size_t k = n; k > 0; --k)
			{
				const std::size_t i = k - 1;
				float value = m_Increment[m_Data.blockEdge[first + i]];
				for (std::size_t j = i + 1; j < std::min(n, i + band + 1); ++j)
				{
					value -= m_Data.lower[offset + j * stride + j - i] * m_Increment[m_Data.blockEdge[first + j]];
				}
				m_Increment[m_Data.blockEdge[first + i]] = value / m_Data.lower[offset + i * stride];
			}
		}
		for (std::size_t i = 0; i < m_Data.electric.size(); ++i)
		{
			m_Data.electric[i] += dt * m_Increment[i];
		}
		for (std::size_t i = 0; i < m_Data.sourceEdge.size(); ++i)
		{
			m_Data.electric[m_Data.sourceEdge[i]] += amplitude * m_Data.sourceWeight[i];
		}
		for (std::size_t i = 0; i < m_Data.magnetic.size(); ++i)
		{
			float value = 0;
			for (std::size_t k = m_Data.curlStart[i]; k < m_Data.curlStart[i + 1]; ++k)
			{
				value += m_Data.curlWeight[k] * m_Data.electric[m_Data.curlColumn[k]];
			}
			m_Data.magnetic[i] -= dt * value;
		}
	}

	void ReadFields(CartesianSubgrid &grid) override
	{
		grid.ImportDeviceFields(m_Data.electric.data(), m_Data.electric.size(),
								m_Data.magnetic.data(), m_Data.magnetic.size());
	}

	std::size_t StorageBytes() const override
	{
		std::size_t floats = m_Increment.capacity();
		for (const auto *values : {&m_Data.electric, &m_Data.magnetic, &m_Data.curlWeight,
			 &m_Data.transposeWeight, &m_Data.lower, &m_Data.sourceWeight})
		{
			floats += values->capacity();
		}
		std::size_t integers = 0;
		for (const auto *values : {&m_Data.curlStart, &m_Data.curlColumn, &m_Data.transposeStart,
			 &m_Data.transposeColumn, &m_Data.blockStart, &m_Data.blockEdge,
			 &m_Data.bandwidth, &m_Data.lowerStart, &m_Data.sourceEdge})
		{
			integers += values->capacity();
		}
		return floats * sizeof(float) + integers * sizeof(std::uint32_t);
	}

  private:
	CartesianSubgridDeviceData m_Data;
	std::vector<float> m_Increment;
};
} // namespace

static std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridReference(const CartesianSubgrid &grid)
{
	return std::make_unique<Reference_CartesianSubgrid>(grid);
}

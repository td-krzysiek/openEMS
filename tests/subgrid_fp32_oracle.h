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
			const float value = Dot(m_Data.transposeStart[i], m_Data.transposeStart[i + 1],
				m_Data.transposeColumn, m_Data.transposeWeight, m_Data.transposeLow, m_Data.magnetic);
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
			m_Data.electric[i] = Update(m_Data.electric[i], timestep, m_Increment[i]);
		}
		for (std::size_t i = 0; i < m_Data.sourceEdge.size(); ++i)
		{
			m_Data.electric[m_Data.sourceEdge[i]] += amplitude * m_Data.sourceWeight[i];
		}
		for (std::size_t i = 0; i < m_Data.magnetic.size(); ++i)
		{
			const float value = Dot(m_Data.curlStart[i], m_Data.curlStart[i + 1],
				m_Data.curlColumn, m_Data.curlWeight, m_Data.curlLow, m_Data.electric);
			m_Data.magnetic[i] = Update(m_Data.magnetic[i], -timestep, value);
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
			 &m_Data.transposeWeight, &m_Data.curlLow, &m_Data.transposeLow, &m_Data.lower, &m_Data.sourceWeight})
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
	// FP32 arithmetic oracle, not a second double-state mathematical reference.
	// The double CPU core is compared independently by the calling test.
	static float Update(float previous, double timestep, float increment)
	{
		const float high = float(timestep), low = float(timestep - double(high));
		const float product = high * increment;
		const float error = std::fma(high, increment, -product) + low * increment;
		float sum = previous, correction = 0;
		for (const float term : {product, error})
		{
			const float next = sum + term;
			correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
			sum = next;
		}
		return sum + correction;
	}
	// Test-only arithmetic oracle; do not call shared device operations here.
	static float Dot(std::size_t begin, std::size_t end, const std::vector<std::uint32_t> &column,
		const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &field)
	{
		float sum = 0, correction = 0;
		for (std::size_t k = begin; k < end; ++k)
		{
			const float input = field[column[k]], product = high[k] * input;
			const float residual = std::fma(high[k], input, -product) + low[k] * input;
			for (const float term : {product, residual})
			{
				const float next = sum + term;
				correction += std::abs(sum) >= std::abs(term) ? (sum - next) + term : (term - next) + sum;
				sum = next;
			}
		}
		return sum + correction;
	}
	CartesianSubgridDeviceData m_Data;
	std::vector<float> m_Increment;
};
} // namespace

static std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridReference(const CartesianSubgrid &grid)
{
	return std::make_unique<Reference_CartesianSubgrid>(grid);
}

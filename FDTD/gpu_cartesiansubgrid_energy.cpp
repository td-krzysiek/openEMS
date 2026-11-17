/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid_energy.h"
#include "gpu_backend.h"
#include "kernels/subgrid_energy.h"
#include "kernels/subgrid_conductivity.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace
{
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}

class Energy : public GPU_CartesianSubgridEnergy
{
  public:
	Energy(GPU_Backend &backend, GPU_CartesianSubgridFields fields,
		const std::vector<CartesianSubgridEnergyCoupling> &couplings)
		: m_Backend(backend)
	{
		const auto limit = std::numeric_limits<std::uint32_t>::max();
		Require(fields.electricCount <= limit && fields.magneticCount <= limit &&
			couplings.size() <= limit - fields.electricCount,
			"CartesianSubgrid: device energy count exceeds indexing");
		Require((!fields.electricCount || fields.electric) && (!fields.magneticCount || fields.magnetic),
			"CartesianSubgrid: device energy requires owned field buffers");
		std::vector<std::uint32_t> first, second;
		std::vector<float> weight, low;
		std::pair<std::uint32_t, std::uint32_t> previous{};
		for (const auto &term : couplings)
		{
			const auto pair = std::make_pair(term.first, term.second);
			Require(term.first < term.second && term.second < fields.electricCount &&
				(first.empty() || previous < pair) && std::isfinite(term.weight) && std::abs(term.weight) < 1,
				"CartesianSubgrid: invalid/noncanonical owned energy coupling");
			const float high = static_cast<float>(term.weight), residual = static_cast<float>(term.weight - double(high));
			Require((term.weight == 0 || std::abs(high) >= std::numeric_limits<float>::min()) &&
				(term.weight == double(high) || std::abs(residual) >= std::numeric_limits<float>::min()),
				"CartesianSubgrid: energy coupling is not representable as normal FP32 pair");
			first.push_back(term.first);
			second.push_back(term.second);
			weight.push_back(high);
			low.push_back(residual);
			previous = pair;
		}
		m_Backend.Activate();
		m_Input = {fields.electric, fields.magnetic};
		m_Counts = {static_cast<std::uint32_t>(fields.electricCount), static_cast<std::uint32_t>(fields.magneticCount)};
		m_First = Upload(first);
		m_Second = Upload(second);
		m_Weight = Upload(weight);
		m_Low = Upload(low);
		m_Pairs = static_cast<std::uint32_t>(couplings.size());
		for (unsigned int field = 0; field < 2; ++field)
		{
			std::size_t count = m_Counts[field];
			while (count)
			{
				count = (count + 511) / 512;
				m_Maximum[field].push_back(Upload(std::vector<float>(count, 0)));
				if (count == 1)
				{
					break;
				}
			}
			count = m_Counts[field] + (field ? 0 : std::size_t(m_Pairs));
			while (count)
			{
				count = (count + 511) / 512;
				m_Sum[field].push_back(Upload(std::vector<float>(2 * count, 0)));
				if (count == 1)
				{
					break;
				}
			}
		}
	}

	std::array<double, 2> Read() override
	{
		m_Backend.Activate();
		for (unsigned int field = 0; field < 2; ++field)
		{
			const float *input = m_Input[field];
			std::size_t count = m_Counts[field];
			for (auto *output : m_Maximum[field])
			{
				Op_SubgridEnergyMax op{input, output, static_cast<std::uint32_t>(count)};
				count = (count + 511) / 512;
				m_Backend.RunPoints(GPU_OP_SUBGRID_ENERGY_MAX, &op, 0, count, false);
				input = output;
			}
		}
		std::array<float, 2> scale{};
		for (unsigned int field = 0; field < 2; ++field)
		{
			if (!m_Maximum[field].empty())
			{
				Download(&scale[field], m_Maximum[field].back(), sizeof(float));
			}
			Require(std::isfinite(scale[field]) && scale[field] >= 0,
				"CartesianSubgrid: nonfinite device field-energy scale");
		}
		for (unsigned int field = 0; field < 2; ++field)
		{
			if (!scale[field])
			{
				continue;
			}
			std::size_t count = m_Counts[field] + (field ? 0 : std::size_t(m_Pairs));
			Op_SubgridEnergyTerms terms{m_Input[field], m_First, m_Second, m_Weight, m_Low,
				m_Sum[field].front(), m_Counts[field], field ? 0 : m_Pairs, scale[field]};
			count = (count + 511) / 512;
			m_Backend.RunPoints(GPU_OP_SUBGRID_ENERGY_TERMS, &terms, 0, count, false);
			for (std::size_t level = 1; level < m_Sum[field].size(); ++level)
			{
				Op_SubgridLossReduce op{m_Sum[field][level - 1], m_Sum[field][level], static_cast<std::uint32_t>(count)};
				count = (count + 511) / 512;
				m_Backend.RunPoints(GPU_OP_SUBGRID_LOSS_REDUCE, &op, 0, count, false);
			}
		}
		std::array<double, 2> result{};
		for (unsigned int field = 0; field < 2; ++field)
		{
			if (scale[field])
			{
				std::array<float, 2> sum{};
				Download(sum.data(), m_Sum[field].back(), sizeof(sum));
				result[field] = (double(sum[0]) - sum[1]) * double(scale[field]) * scale[field];
				Require(std::isfinite(result[field]) && result[field] >= 0,
					"CartesianSubgrid: nonfinite/negative device field energy");
			}
		}
		return result;
	}
	std::size_t StorageBytes() const override { return m_Bytes; }
	std::uint64_t DownloadBytes() const override { return m_DownloadBytes; }

  private:
	template <class T> T *Upload(const std::vector<T> &values)
	{
		const auto bytes = std::max(values.size() * sizeof(T), std::size_t(4));
		void *result = m_Backend.DeviceAlloc(bytes, values.empty() ? nullptr : values.data());
		if (!result)
		{
			throw std::runtime_error("CartesianSubgrid: device energy allocation failed");
		}
		m_Bytes += bytes;
		return static_cast<T *>(result);
	}
	void Download(void *destination, const void *source, std::size_t bytes)
	{
		m_Backend.DeviceDownload(destination, source, bytes);
		m_DownloadBytes += bytes;
	}
	GPU_Backend &m_Backend;
	std::array<const float *, 2> m_Input{};
	std::array<std::uint32_t, 2> m_Counts{};
	std::uint32_t m_Pairs = 0;
	std::uint32_t *m_First = nullptr, *m_Second = nullptr;
	float *m_Weight = nullptr, *m_Low = nullptr;
	std::array<std::vector<float *>, 2> m_Maximum, m_Sum;
	std::size_t m_Bytes = 0;
	std::uint64_t m_DownloadBytes = 0;
};
}

std::unique_ptr<GPU_CartesianSubgridEnergy> GPU_CreateCartesianSubgridEnergy(
	GPU_Backend &backend, GPU_CartesianSubgridFields fields,
	const std::vector<CartesianSubgridEnergyCoupling> &couplings)
{
	return std::make_unique<Energy>(backend, fields, couplings);
}

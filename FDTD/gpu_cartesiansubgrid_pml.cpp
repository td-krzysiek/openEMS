/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid_pml.h"
#include "kernels/subgrid_pml.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
class CompactDevicePML : public GPU_CartesianSubgridPML
{
  public:
	CompactDevicePML(GPU_Backend &backend, const GPU_CartesianSubgridFields &fields,
		const CartesianSubgridPMLData &data)
		: m_Backend(backend)
	{
		if (!fields.electric || !fields.magnetic || !(data.timestep > 0) || !std::isfinite(data.timestep) ||
			fields.fixedTimestep != data.timestep)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid paired device PML fields/timestep");
		}
		m_Backend.Activate();
		Init(m_Electric, data.electric, fields.magnetic, fields.magneticCount,
			fields.electric, fields.electricCount, data.timestep);
		Init(m_Magnetic, data.magnetic, fields.electric, fields.electricCount,
			fields.magnetic, fields.magneticCount, -data.timestep);
	}

	void Hook(bool currents, GPU_Hook hook, const Engine_Box &, unsigned int ts) override
	{
		const unsigned int phase = (currents ? 3 : 0) + static_cast<unsigned int>(hook);
		if (phase != m_Phase || (m_Started && ts != m_Timestep) || (!m_Started && ts != 0))
		{
			throw std::invalid_argument("CartesianSubgrid: device PML hook order/timestep mismatch");
		}
		// Another worker may have selected a different physical device since
		// the preceding hook. Fields and stream still belong to this backend.
		m_Backend.Activate();
		if (!m_Started)
		{
			m_Timestep = ts;
			m_Started = true;
		}
		auto &state = currents ? m_Magnetic : m_Electric;
		if (hook == GPU_PRE)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_PML_PRE, &state.pre, 0, state.rows, false);
		}
		else
		{
			state.post.apply = hook == GPU_APPLY ? 1 : 0;
			m_Backend.RunPoints(GPU_OP_SUBGRID_PML_POST, &state.post, 0, state.groups, false);
		}
		m_Phase = (m_Phase + 1) % 6;
		if (m_Phase == 0)
		{
			++m_Timestep;
		}
	}

	double ReadBoundaryWorkOut() override
	{
		if (m_Phase != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: PML diagnostics require a complete timestep");
		}
		m_Backend.Activate();
		double total = 0;
		for (auto *state : {&m_Electric, &m_Magnetic})
		{
			if (!state->work.empty())
			{
				m_Backend.DeviceDownload(state->work.data(), state->post.work, state->work.size() * sizeof(float));
			}
			for (std::size_t i = 0; i < state->work.size(); i += 2)
			{
				total += double(state->work[i]) - state->work[i + 1];
			}
		}
		if (!std::isfinite(total))
		{
			throw std::runtime_error("CartesianSubgrid: nonfinite device PML work");
		}
		return total;
	}

	std::size_t StorageBytes() const override
	{
		return m_StorageBytes;
	}

  private:
	struct State
	{
		Op_SubgridPMLPre pre{};
		Op_SubgridPMLPost post{};
		std::vector<float> work;
		std::size_t rows = 0, groups = 0;
	};

	void Init(State &state, const CartesianSubgridPMLData::Rows &data,
		const float *input, std::size_t sourceCount, float *field, std::size_t targetCount, double timestep)
	{
		const std::size_t rows = data.target.size();
		if (data.start.size() != rows + 1 || data.decay.size() != rows || data.change.size() != rows ||
			data.weight.size() != data.source.size() || data.low.size() != data.source.size() ||
			data.start.front() != 0 || data.start.back() != data.source.size() ||
			data.group.empty() || data.group.front() != 0 || data.group.back() != rows ||
			data.group.size() - 1 > std::numeric_limits<std::uint32_t>::max() / 2)
		{
			throw std::invalid_argument("CartesianSubgrid: malformed device PML tables");
		}
		for (std::size_t row = 0; row < rows; ++row)
		{
			if (data.target[row] >= targetCount || data.start[row] >= data.start[row + 1] ||
				!(data.decay[row] > 0 && data.decay[row] <= 1) ||
				!(data.change[row] >= -1 && data.change[row] <= 0) ||
				data.decay[row] < std::numeric_limits<float>::min() ||
				(data.change[row] != 0 && std::abs(data.change[row]) < std::numeric_limits<float>::min()))
			{
				throw std::invalid_argument("CartesianSubgrid: invalid device PML row/rate");
			}
		}
		for (std::size_t k = 0; k < data.source.size(); ++k)
		{
			if (data.source[k] >= sourceCount || !std::isfinite(data.weight[k]) || !std::isfinite(data.low[k]) ||
				(data.weight[k] != 0 && std::abs(data.weight[k]) < std::numeric_limits<float>::min()) ||
				(data.low[k] != 0 && std::abs(data.low[k]) < std::numeric_limits<float>::min()))
			{
				throw std::invalid_argument("CartesianSubgrid: invalid device PML derivative");
			}
		}
		for (std::size_t group = 0; group + 1 < data.group.size(); ++group)
		{
			const std::size_t first = data.group[group], stop = data.group[group + 1];
			if (first >= stop || stop > rows || stop - first > 3 ||
				(first && data.target[first - 1] >= data.target[first]))
			{
				throw std::invalid_argument("CartesianSubgrid: invalid PML target grouping");
			}
			for (std::size_t row = first; row < stop; ++row)
			{
				if (data.target[row] != data.target[first])
				{
					throw std::invalid_argument("CartesianSubgrid: PML group has different targets");
				}
			}
		}
		const float high = float(timestep), low = float(timestep - double(high));
		if (!std::isfinite(high) || std::abs(high) < std::numeric_limits<float>::min() ||
			(timestep != double(high) && std::abs(low) < std::numeric_limits<float>::min()))
		{
			throw std::invalid_argument("CartesianSubgrid: unrepresentable device PML timestep");
		}
		state.rows = rows;
		state.groups = data.group.size() - 1;
		state.work.assign(2 * state.groups, 0);
		auto *target = Upload(data.target);
		auto *correction = Upload(std::vector<float>(rows, 0));
		auto *previous = Upload(std::vector<float>(rows, 0));
		state.pre = {Upload(data.start), Upload(data.source), target, Upload(data.weight), Upload(data.low),
			Upload(data.decay), Upload(data.change), input, field, Upload(std::vector<float>(rows, 0)),
			correction, previous, high, low};
		state.post = {Upload(data.group), target, correction, previous, field,
			Upload(state.work), Upload(std::vector<float>(state.groups, 0)), 0};
	}

	template <class T> T *Upload(const std::vector<T> &values)
	{
		const std::size_t bytes = std::max(values.size() * sizeof(T), std::size_t(4));
		void *buffer = m_Backend.DeviceAlloc(bytes, values.empty() ? nullptr : values.data());
		if (!buffer)
		{
			throw std::runtime_error("CartesianSubgrid: device PML allocation failed");
		}
		m_StorageBytes += bytes;
		return static_cast<T *>(buffer);
	}

	GPU_Backend &m_Backend;
	State m_Electric, m_Magnetic;
	std::size_t m_StorageBytes = 0;
	unsigned int m_Phase = 0, m_Timestep = 0;
	bool m_Started = false;
};
} // namespace

std::unique_ptr<GPU_CartesianSubgridPML> GPU_CreateCartesianSubgridPML(
	GPU_Backend &backend, const GPU_CartesianSubgridFields &fields, const CartesianSubgridPMLData &data)
{
	return std::make_unique<CompactDevicePML>(backend, fields, data);
}

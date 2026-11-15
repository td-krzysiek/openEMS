/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid.h"
#include "gpu_backend.h"
#include "cartesian_subgrid.h"
#include "kernels/subgrid_update.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace
{
//! Existing backend allocation/point operations; fields stay on that backend.
class CartesianSubgridExecution : public GPU_CartesianSubgrid
{
  public:
	CartesianSubgridExecution(GPU_Backend &backend, const CartesianSubgrid &grid)
		: m_Backend(backend)
	{
		const auto data = grid.DeviceData();
		m_StableStep = data.stableTimestep;
		m_Electric.resize(data.electric.size());
		m_Magnetic.resize(data.magnetic.size());
		m_BlockCount = data.bandwidth.size();
		m_SourceCount = data.sourceEdge.size();
		auto *electric = Upload(data.electric), *magnetic = Upload(data.magnetic);
		auto *increment = Upload(std::vector<float>(data.electric.size(), 0));
		m_Transpose = {Upload(data.transposeStart), Upload(data.transposeColumn),
			Upload(data.transposeWeight), magnetic, increment, 1, 0};
		m_Curl = {Upload(data.curlStart), Upload(data.curlColumn),
			Upload(data.curlWeight), electric, magnetic, 0, 1};
		m_Mass = {Upload(data.blockStart), Upload(data.blockEdge), Upload(data.bandwidth),
			Upload(data.lowerStart), Upload(data.lower), increment};
		m_Update = {electric, increment, 0};
		m_Source = {electric, Upload(data.sourceEdge), Upload(data.sourceWeight), 0};
	}

	~CartesianSubgridExecution() override
	{
		// Buffers are owned by the creating backend, whose destructor waits for
		// queued work before releasing them. Do not throw or add a wait here.
	}

	void Step(double timestep, double sourceAmplitude) override
	{
		const float dt = static_cast<float>(timestep), amplitude = static_cast<float>(sourceAmplitude);
		if (!(timestep > 0) || timestep > m_StableStep || !std::isfinite(dt) ||
			dt < std::numeric_limits<float>::min() || !std::isfinite(amplitude))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid device step or drive");
		}
		m_Update.timestep = dt;
		m_Curl.scale = -dt;
		m_Source.amplitude = amplitude;
		m_Backend.RunPoints(GPU_OP_SUBGRID_CURL, &m_Transpose, 0, m_Electric.size(), false);
		m_Backend.RunPoints(GPU_OP_SUBGRID_MASS, &m_Mass, 0, m_BlockCount, false);
		m_Backend.RunPoints(GPU_OP_SUBGRID_ELECTRIC, &m_Update, 0, m_Electric.size(), false);
		if (m_SourceCount && amplitude != 0)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_SOURCE, &m_Source, 0, m_SourceCount, true);
		}
		m_Backend.RunPoints(GPU_OP_SUBGRID_CURL, &m_Curl, 0, m_Magnetic.size(), false);
	}

	void ReadFields(CartesianSubgrid &grid) override
	{
		m_Backend.DeviceDownload(m_Electric.data(), m_Update.electric, m_Electric.size() * sizeof(float));
		m_Backend.DeviceDownload(m_Magnetic.data(), m_Curl.output, m_Magnetic.size() * sizeof(float));
		grid.ImportDeviceFields(m_Electric.data(), m_Electric.size(), m_Magnetic.data(), m_Magnetic.size());
	}

	std::size_t StorageBytes() const override
	{
		return m_StorageBytes;
	}

  private:
	template <class T> T *Upload(const std::vector<T> &values)
	{
		const std::size_t bytes = std::max(values.size() * sizeof(T), std::size_t(4));
		void *buffer = m_Backend.DeviceAlloc(bytes, values.empty() ? nullptr : values.data());
		if (!buffer)
		{
			throw std::runtime_error("CartesianSubgrid: device allocation failed");
		}
		m_StorageBytes += bytes;
		return static_cast<T *>(buffer);
	}

	GPU_Backend &m_Backend;
	std::size_t m_StorageBytes = 0, m_BlockCount = 0, m_SourceCount = 0;
	double m_StableStep = 0;
	std::vector<float> m_Electric, m_Magnetic;
	Op_SubgridCurl m_Transpose{}, m_Curl{};
	Op_SubgridMass m_Mass{};
	Op_SubgridElectric m_Update{};
	Op_SubgridSource m_Source{};
};
} // namespace

std::unique_ptr<GPU_CartesianSubgrid> GPU_Backend::CreateCartesianSubgrid(const CartesianSubgrid &grid)
{
	return std::make_unique<CartesianSubgridExecution>(*this, grid);
}

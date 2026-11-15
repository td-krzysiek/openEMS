/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid.h"
#include "gpu_backend.h"
#include "cartesian_subgrid.h"
#include "kernels/subgrid_update.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <iostream>
#include <sstream>
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
		if (data.factorTails.terms)
		{
			std::ostringstream message;
			message << std::scientific << "CartesianSubgrid FP32 factor FTZ tails: " << data.factorTails.terms
					<< "; dropped L1=" << data.factorTails.l1
					<< "; maximum relative factor bound=" << data.factorTails.relativeBound;
			std::cout << message.str() << std::endl;
		}
		if (data.stateTails.terms)
		{
			std::ostringstream message;
			message << std::scientific << "CartesianSubgrid FP32 initial state FTZ tails: " << data.stateTails.terms
				<< "; dropped L1=" << data.stateTails.l1
				<< "; maximum relative state bound=" << data.stateTails.relativeBound;
			std::cout << message.str() << std::endl;
		}
		m_StableStep = data.stableTimestep;
		m_PortStep = data.portTimestep;
		m_PortCount = data.rootConductance.size();
		m_Electric.resize(data.electric.size());
		m_Magnetic.resize(data.magnetic.size());
		m_BlockCount = data.bandwidth.size();
		m_SourceCount = data.sourceEdge.size();
		auto *electric = Upload(data.electric), *magnetic = Upload(data.magnetic);
		auto *increment = Upload(std::vector<float>(data.electric.size(), 0));
		m_Transpose = {Upload(data.transposeStart), Upload(data.transposeColumn),
			Upload(data.transposeWeight), Upload(data.transposeLow), magnetic, increment, 1, 0, 0};
		m_Curl = {Upload(data.curlStart), Upload(data.curlColumn),
			Upload(data.curlWeight), Upload(data.curlLow), electric, magnetic, 0, 0, 1};
		m_Mass = {Upload(data.blockStart), Upload(data.blockEdge), Upload(data.bandwidth),
			Upload(data.lowerStart), Upload(data.lower), increment};
		m_Update = {electric, increment, 0, 0};
		m_LedgerDevice = Upload(data.ledger);
		auto *sourceEdge = Upload(data.sourceEdge);
		auto *sourceCurl = Upload(std::vector<float>(data.sourceEdge.size(), 0));
		m_SourceCurl = {increment, sourceEdge, sourceCurl};
		m_Source = {electric, sourceEdge, Upload(data.sourceWeight), Upload(data.sourceOffStart),
			Upload(data.sourceOffEdge), Upload(data.sourceOffWeight), sourceCurl, m_LedgerDevice,
			0, 0, static_cast<std::uint32_t>(m_PortCount != 0)};
		auto *portRHS = Upload(std::vector<float>(m_PortCount, 0));
		m_PortSolve = {Upload(data.portStart), Upload(data.portEdge), Upload(data.portWeight),
			Upload(data.rootConductance), Upload(data.portCurrent), Upload(data.portLower),
			electric, increment, portRHS, m_LedgerDevice, 0, 0, static_cast<std::uint32_t>(m_PortCount)};
		m_PortApply = {electric, Upload(data.responseStart), Upload(data.responseEdge),
			Upload(data.responseWeight), portRHS, 0};
	}

	~CartesianSubgridExecution() override
	{
		// Buffers are owned by the creating backend, whose destructor waits for
		// queued work before releasing them. Do not throw or add a wait here.
	}

	void Step(double timestep, double sourceAmplitude) override
	{
		UpdateElectric(timestep, sourceAmplitude);
		UpdateMagnetic(timestep);
	}

	void RequireFixedTimestep(double timestep) override
	{
		if (m_Started || !(timestep > 0) || timestep > m_StableStep ||
			(m_PortCount && timestep != m_PortStep) || (m_FixedStep != 0 && timestep != m_FixedStep))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid fixed device recurrence timestep");
		}
		m_FixedStep = timestep;
	}

	GPU_CartesianSubgridFields Fields() const override
	{
		return {m_Update.electric, m_Curl.output, m_Electric.size(), m_Magnetic.size(),
			m_FixedStep != 0 ? m_FixedStep : m_PortStep};
	}

	void UpdateElectric(double timestep, double sourceAmplitude) override
	{
		if (m_PendingStep != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: electric update requires a complete prior step");
		}
		const float dt = static_cast<float>(timestep), amplitude = static_cast<float>(sourceAmplitude);
		const double residual = timestep - double(dt);
		const float dtLow = static_cast<float>(residual);
		if (!(timestep > 0) || timestep > m_StableStep || !std::isfinite(dt) ||
			dt < std::numeric_limits<float>::min() || !std::isfinite(amplitude) ||
			(residual != 0 && std::abs(dtLow) < std::numeric_limits<float>::min()) ||
			((m_PortCount && timestep != m_PortStep) || (m_FixedStep != 0 && timestep != m_FixedStep)))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid device step or drive");
		}
		m_Update.timestep = dt;
		m_Update.timestepLow = dtLow;
		m_Curl.scale = -dt;
		m_Curl.scaleLow = -m_Update.timestepLow;
		m_Source.amplitude = amplitude;
		m_Source.timestep = dt;
		m_PortSolve.timestep = dt;
		m_PortSolve.amplitude = amplitude;
		m_PortApply.timestep = dt;
		m_Backend.RunPoints(GPU_OP_SUBGRID_CURL, &m_Transpose, 0, m_Electric.size(), false);
		if (m_PortCount && m_SourceCount)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_SOURCE_CURL, &m_SourceCurl, 0, m_SourceCount, false);
		}
		m_Backend.RunPoints(GPU_OP_SUBGRID_MASS, &m_Mass, 0, m_BlockCount, false);
		if (m_PortCount)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_PORT_SOLVE, &m_PortSolve, 0, 1, false);
		}
		m_Backend.RunPoints(GPU_OP_SUBGRID_ELECTRIC, &m_Update, 0, m_Electric.size(), false);
		if (m_PortCount)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_PORT_APPLY, &m_PortApply, 0, m_PortCount, true);
		}
		if (m_SourceCount && amplitude != 0)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_SOURCE, &m_Source, 0, m_SourceCount, true);
		}
		m_Started = true;
		m_PendingStep = timestep;
	}

	void UpdateMagnetic(double timestep) override
	{
		if (m_PendingStep == 0 || timestep != m_PendingStep)
		{
			throw std::invalid_argument("CartesianSubgrid: magnetic update requires its paired electric step");
		}
		m_Backend.RunPoints(GPU_OP_SUBGRID_CURL, &m_Curl, 0, m_Magnetic.size(), false);
		m_PendingStep = 0;
	}

	void ReadFields(CartesianSubgrid &grid) override
	{
		if (m_PendingStep != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: device observation requires a complete timestep");
		}
		m_Backend.DeviceDownload(m_Electric.data(), m_Update.electric, m_Electric.size() * sizeof(float));
		m_Backend.DeviceDownload(m_Magnetic.data(), m_Curl.output, m_Magnetic.size() * sizeof(float));
		grid.ImportDeviceFields(m_Electric.data(), m_Electric.size(), m_Magnetic.data(), m_Magnetic.size());
		if (m_PortCount)
		{
			m_Backend.DeviceDownload(m_Ledger.data(), m_LedgerDevice, sizeof(m_Ledger));
			grid.ImportDeviceLedger(m_Ledger[0], m_Ledger[2]);
		}
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
	std::size_t m_StorageBytes = 0, m_BlockCount = 0, m_SourceCount = 0, m_PortCount = 0;
	double m_StableStep = 0, m_PortStep = 0, m_FixedStep = 0, m_PendingStep = 0;
	bool m_Started = false;
	std::vector<float> m_Electric, m_Magnetic;
	Op_SubgridCurl m_Transpose{}, m_Curl{};
	Op_SubgridMass m_Mass{};
	Op_SubgridElectric m_Update{};
	Op_SubgridSource m_Source{};
	Op_SubgridSourceCurl m_SourceCurl{};
	Op_SubgridPortSolve m_PortSolve{};
	Op_SubgridPortApply m_PortApply{};
	float *m_LedgerDevice = nullptr;
	std::array<float, 4> m_Ledger{};
};
} // namespace

std::unique_ptr<GPU_CartesianSubgrid> GPU_Backend::CreateCartesianSubgrid(const CartesianSubgrid &grid)
{
	return std::make_unique<CartesianSubgridExecution>(*this, grid);
}

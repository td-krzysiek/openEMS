/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid_conductivity.h"
#include "gpu_backend.h"
#include "subgrid_device_curl.h"
#include "kernels/subgrid_conductivity.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace
{
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}
float Normal(double value)
{
	const float result = static_cast<float>(value);
	Require(std::isfinite(value) && std::isfinite(result) &&
		(value == 0 || std::abs(result) >= std::numeric_limits<float>::min()),
		"CartesianSubgrid: unrepresentable device conductivity coefficient");
	return result;
}

class CompactDeviceConductivity : public GPU_CartesianSubgridRecurrence
{
  public:
	explicit CompactDeviceConductivity(CartesianSubgridConductivityData data) : m_Data(std::move(data)) {}
	void Bind(GPU_Backend &backend, const GPU_CartesianSubgridFields &fields, float *rhs) override
	{
		Require(!m_Backend && fields.electric && fields.magnetic && rhs &&
			m_Data.timestep > 0 && std::isfinite(m_Data.timestep) &&
			m_Data.magneticRate.size() == fields.magneticCount &&
			m_Data.conductance.start.size() == fields.electricCount + 1,
			"CartesianSubgrid: conductivity device bind mismatch");
		m_Backend = &backend;
		m_Timestep = m_Data.timestep;
		backend.Activate();
		const auto &g = m_Data.conductance;
		const CartesianSubgridDeviceCurl conductance(g.start, g.edge, g.weight, g.low, fields.electricCount);
		m_Conductance = {Upload(conductance.start), Upload(conductance.offset), Upload(conductance.base),
			Upload(conductance.pattern), Upload(conductance.weight), Upload(conductance.residual),
			fields.electric, rhs, -1, 0, 1};
		const auto &t = m_Data.trace;
		// Validate CSR bounds/coefficients with the same portable row contract.
		const CartesianSubgridDeviceCurl validated(t.start, t.edge, t.weight, t.low, fields.electricCount, true);
		(void)validated;
		m_ElectricCount = fields.electricCount;
		m_MagneticCount = fields.magneticCount;
		m_TraceCount = t.start.size() - 1;
		Require(m_TraceCount <= std::numeric_limits<std::uint32_t>::max() / 2 &&
			m_MagneticCount <= std::numeric_limits<std::uint32_t>::max() / 2,
			"CartesianSubgrid: conductivity ledger exceeds device indexing");
		if (m_TraceCount)
		{
			std::vector<std::uint32_t> edges(fields.electricCount);
			std::iota(edges.begin(), edges.end(), 0);
			m_Previous = {fields.electric, Upload(edges), Upload(std::vector<float>(fields.electricCount, 0))};
			m_ElectricLoss = {Upload(t.start), Upload(t.edge), Upload(t.weight), Upload(t.low),
				m_Previous.curl, fields.electric, Upload(std::vector<float>(2 * m_TraceCount, 0)), Normal(m_Timestep)};
			InitReduction(m_ElectricReduction, m_TraceCount);
		}
		std::vector<float> decay, decayLow, change, changeLow, rateStep;
		for (double rate : m_Data.magneticRate)
		{
			Require(rate >= 0 && std::isfinite(rate) && std::isfinite(m_Timestep * rate),
				"CartesianSubgrid: invalid device magnetic loss rate");
			m_HasMagneticLoss = m_HasMagneticLoss || rate != 0;
			const double a = .5 * m_Timestep * rate;
			const double d = -m_Timestep * rate / (1 + a), c = -m_Timestep / (1 + a);
			decay.push_back(Normal(d));
			decayLow.push_back(Normal(d - double(decay.back())));
			change.push_back(Normal(c));
			changeLow.push_back(Normal(c - double(change.back())));
			rateStep.push_back(Normal(m_Timestep * rate));
		}
		if (m_HasMagneticLoss)
		{
			m_MagneticLoss = {fields.magnetic, Upload(std::vector<float>(2 * fields.magneticCount, 0)),
				Upload(std::vector<float>(fields.magneticCount, 0)), Upload(decay), Upload(decayLow),
				Upload(change), Upload(changeLow), Upload(rateStep)};
			InitReduction(m_MagneticReduction, fields.magneticCount);
		}
		m_Data = {};
	}
	void PrepareElectric(double timestep) override
	{
		Require(m_Backend && m_Phase == 0 && timestep == m_Timestep,
			"CartesianSubgrid: conductivity electric phase/timestep mismatch");
		m_Phase = 1;
		if (m_TraceCount)
		{
			m_Backend->RunPoints(GPU_OP_SUBGRID_GATHER, &m_Previous, 0, m_ElectricCount, false);
		}
		m_Backend->RunPoints(GPU_OP_SUBGRID_CURL, &m_Conductance, 0, m_ElectricCount, false);
	}
	void FinishElectric(double timestep) override
	{
		Require(m_Phase == 1 && timestep == m_Timestep, "CartesianSubgrid: conductivity implicit phase mismatch");
		m_Phase = 2;
		if (m_TraceCount)
		{
			m_Backend->RunPoints(GPU_OP_SUBGRID_ELECTRIC_LOSS, &m_ElectricLoss, 0, m_TraceCount, false);
		}
	}
	void UpdateMagnetic(const Op_SubgridCurl &curl, double timestep) override
	{
		Require(m_Phase == 2 && timestep == m_Timestep, "CartesianSubgrid: conductivity magnetic phase mismatch");
		m_Phase = 3;
		if (m_HasMagneticLoss)
		{
			auto temporary = curl;
			temporary.output = const_cast<float *>(m_MagneticLoss.curl);
			temporary.additive = 0;
			m_Backend->RunPoints(GPU_OP_SUBGRID_CURL, &temporary, 0, m_MagneticCount, false);
			m_Backend->RunPoints(GPU_OP_SUBGRID_MAGNETIC_LOSS, &m_MagneticLoss, 0, m_MagneticCount, false);
		}
		else
		{
			m_Backend->RunPoints(GPU_OP_SUBGRID_CURL, &curl, 0, m_MagneticCount, false);
		}
		m_Phase = 0;
	}
	void CheckObservation() const override
	{
		Require(m_Backend && m_Phase == 0, "CartesianSubgrid: conductivity observation requires a complete step");
	}
	std::array<double, 2> ReadLedger() override
	{
		CheckObservation();
		return {m_TraceCount ? Reduce(m_ElectricLoss.loss, m_ElectricReduction, m_TraceCount) : 0,
			m_HasMagneticLoss ? Reduce(m_MagneticLoss.loss, m_MagneticReduction, m_MagneticCount) : 0};
	}
	std::size_t StorageBytes() const override { return m_StorageBytes; }

  private:
	void InitReduction(std::vector<float *> &buffers, std::size_t count)
	{
		Require(count <= std::numeric_limits<std::uint32_t>::max() / 2,
			"CartesianSubgrid: conductivity ledger exceeds device indexing");
		do
		{
			count = (count + 511) / 512;
			buffers.push_back(Upload(std::vector<float>(2 * count, 0)));
		} while (count > 1);
	}
	double Reduce(const float *input, const std::vector<float *> &buffers, std::size_t count)
	{
		for (float *output : buffers)
		{
			Op_SubgridLossReduce op{input, output, static_cast<std::uint32_t>(count)};
			count = (count + 511) / 512;
			m_Backend->RunPoints(GPU_OP_SUBGRID_LOSS_REDUCE, &op, 0, count, false);
			input = output;
		}
		std::array<float, 2> result{};
		m_Backend->DeviceDownload(result.data(), input, sizeof(result));
		const double value = double(result[0]) - result[1];
		Require(value >= 0 && std::isfinite(value), "CartesianSubgrid: nonfinite/negative device conductivity loss");
		return value;
	}
	template <class T> T *Upload(const std::vector<T> &values)
	{
		const auto bytes = std::max(values.size() * sizeof(T), std::size_t(4));
		void *buffer = m_Backend->DeviceAlloc(bytes, values.empty() ? nullptr : values.data());
		if (!buffer)
		{
			throw std::runtime_error("CartesianSubgrid: device conductivity allocation failed");
		}
		m_StorageBytes += bytes;
		return static_cast<T *>(buffer);
	}
	CartesianSubgridConductivityData m_Data;
	GPU_Backend *m_Backend = nullptr;
	Op_SubgridCurl m_Conductance{};
	Op_SubgridGather m_Previous{};
	Op_SubgridElectricLoss m_ElectricLoss{};
	Op_SubgridMagneticLoss m_MagneticLoss{};
	std::vector<float *> m_ElectricReduction, m_MagneticReduction;
	std::size_t m_StorageBytes = 0, m_ElectricCount = 0, m_MagneticCount = 0, m_TraceCount = 0;
	double m_Timestep = 0;
	unsigned int m_Phase = 0;
	bool m_HasMagneticLoss = false;
};
}
std::unique_ptr<GPU_CartesianSubgridRecurrence> GPU_CreateCartesianSubgridConductivity(CartesianSubgridConductivityData data)
{
	return std::make_unique<CompactDeviceConductivity>(std::move(data));
}

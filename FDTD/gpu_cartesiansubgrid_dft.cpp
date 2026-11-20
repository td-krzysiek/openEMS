/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid_dft.h"
#include "gpu_backend.h"
#include "kernels/subgrid_dft.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}

class DFT : public GPU_CartesianSubgridDFT
{
  public:
	DFT(GPU_Backend &backend, const float *field, std::size_t owned,
		const std::vector<std::uint32_t> &indices, std::size_t frequencies)
		: m_Backend(backend), m_Bytes(CartesianSubgridDFTBytes(indices.size(), frequencies)), m_Frequencies(frequencies)
	{
		Require((indices.empty() || field) && std::all_of(indices.begin(), indices.end(),
			[&](std::uint32_t index) { return index < owned; }),
			"CartesianSubgrid: DFT requires selected owned states");
		m_Op.field = field;
		m_Op.count = static_cast<std::uint32_t>(indices.size());
		if (!indices.empty())
		{
			m_Backend.Activate();
			m_Op.indices = Upload(indices);
			m_Op.sums = Upload(std::vector<float>(4 * indices.size() * frequencies, 0));
		}
	}
	void Accumulate(const std::vector<std::complex<float>> &weights) override
	{
		CartesianSubgridDFTWeights(weights, m_Frequencies);
		if (!m_Op.count)
		{
			return;
		}
		m_Backend.Activate();
		for (std::size_t k = 0; k < weights.size(); ++k)
		{
			auto op = m_Op;
			op.frequency = static_cast<std::uint32_t>(k);
			op.real = weights[k].real();
			op.imag = weights[k].imag();
			// Inline launch arguments avoid a blocking existing-buffer upload.
			m_Backend.RunPoints(GPU_OP_SUBGRID_DFT, &op, 0, m_Op.count, false);
		}
	}
	void Read(std::size_t frequency, std::vector<std::complex<double>> &sums) override
	{
		Require(frequency < m_Frequencies, "CartesianSubgrid: DFT frequency outside configured group");
		std::vector<float> values(4 * std::size_t(m_Op.count));
		if (!values.empty())
		{
			m_Backend.Activate();
			m_Backend.DeviceDownload(values.data(), m_Op.sums + frequency * values.size(), values.size() * sizeof(float));
			m_DownloadBytes += values.size() * sizeof(float);
		}
		sums.resize(m_Op.count);
		for (std::size_t k = 0; k < sums.size(); ++k)
		{
			const double real = double(values[4 * k]) + values[4 * k + 1];
			const double imag = double(values[4 * k + 2]) + values[4 * k + 3];
			Require(std::isfinite(real) && std::isfinite(imag), "CartesianSubgrid: nonfinite DFT state sum");
			sums[k] = {real, imag};
		}
	}
	std::size_t StorageBytes() const override { return m_Bytes; }
	std::uint64_t DownloadBytes() const override { return m_DownloadBytes; }

  private:
	template <class T> T *Upload(const std::vector<T> &values)
	{
		void *buffer = m_Backend.DeviceAlloc(values.size() * sizeof(T), values.data());
		if (!buffer)
		{
			throw std::runtime_error("CartesianSubgrid: DFT allocation failed");
		}
		return static_cast<T *>(buffer);
	}
	GPU_Backend &m_Backend;
	Op_SubgridDFT m_Op{};
	std::size_t m_Bytes, m_Frequencies;
	std::uint64_t m_DownloadBytes = 0;
};
}

std::size_t CartesianSubgridDFTBytes(std::size_t count, std::size_t frequencies)
{
	const auto indexLimit = std::numeric_limits<std::uint32_t>::max();
	const auto limit = std::numeric_limits<std::size_t>::max();
	Require(count <= indexLimit && frequencies > 0 && frequencies <= indexLimit,
		"CartesianSubgrid: invalid DFT selection/frequency count");
	if (!count)
	{
		return 0;
	}
	Require(count <= limit / 16 && frequencies <= (limit - 4 * count) / (16 * count),
		"CartesianSubgrid: DFT buffer size overflow");
	return 4 * count + frequencies * (16 * count);
}
void CartesianSubgridDFTWeights(const std::vector<std::complex<float>> &weights, std::size_t frequencies)
{
	Require(weights.size() == frequencies && std::all_of(weights.begin(), weights.end(),
		[](const std::complex<float> &weight)
		{
			return std::isfinite(weight.real()) && std::isfinite(weight.imag());
		}), "CartesianSubgrid: DFT requires one finite complex weight per frequency");
}
std::unique_ptr<GPU_CartesianSubgridDFT> GPU_CreateCartesianSubgridDFT(
	GPU_Backend &backend, const float *field, std::size_t owned,
	const std::vector<std::uint32_t> &indices, std::size_t frequencies)
{
	return std::make_unique<DFT>(backend, field, owned, indices, frequencies);
}

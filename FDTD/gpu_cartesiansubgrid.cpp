/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "gpu_cartesiansubgrid.h"
#include "gpu_cartesiansubgrid_conductivity.h"
#include "gpu_cartesiansubgrid_energy.h"
#include "gpu_cartesiansubgrid_dft.h"
#include "subgrid_device_curl.h"
#include "subgrid_device_geometry.h"
#include "gpu_backend.h"
#include "cartesian_subgrid.h"
#include "cartesian_subgrid_domains.h"
#include "subgrid_interface_limits.h"
#include "kernels/subgrid_update.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{
//! The expanded contract remains available for independent controls and domains.
CartesianSubgridDeviceOperators PackExpanded(const CartesianSubgridDeviceData &data, bool domain)
{
	CartesianSubgridDeviceOperators result;
	CartesianSubgridDeviceCurl transpose(data.transposeStart, data.transposeColumn,
		data.transposeWeight, data.transposeLow, data.magnetic.size(), domain);
	CartesianSubgridDeviceCurl curl(data.curlStart, data.curlColumn,
		data.curlWeight, data.curlLow, data.electric.size(), domain);
	const auto spatialRows = domain ? 0 : CartesianSubgridDeviceGeometry::Validate(
		data.curlGeometry, data.curlMetrics, data.electric.size(), data.magnetic.size());
	if (spatialRows)
	{
		CartesianSubgridDeviceGeometry::ValidateRows(data);
		std::vector<std::uint32_t> start{0}, columns;
		std::vector<float> high, low;
		for (std::size_t row = 0; row < data.magnetic.size(); ++row)
		{
			if (!subgrid_geometry_regular(data.curlGeometry.data(), std::uint32_t(row)))
			{
				for (std::size_t k = data.curlStart[row]; k < data.curlStart[row + 1]; ++k)
				{
					columns.push_back(data.curlColumn[k]);
					high.push_back(data.curlWeight[k]);
					low.push_back(data.curlLow[k]);
				}
				start.push_back(static_cast<std::uint32_t>(columns.size()));
			}
		}
		CartesianSubgridDeviceCurl packedCurl(start, columns, high, low, data.electric.size());
		high = data.transposeWeight;
		low = data.transposeLow;
		for (std::size_t edge = 0; edge < data.electric.size(); ++edge)
		{
			for (std::size_t k = data.transposeStart[edge]; k < data.transposeStart[edge + 1]; ++k)
			{
				const auto face = data.transposeColumn[k];
				if (!subgrid_geometry_regular(data.curlGeometry.data(), face))
				{
					continue;
				}
				// A regular transpose term stores an exact 0..3 slot, not an
				// approximate coefficient. Both gathers read the same pair.
				unsigned int term = 0;
				while (term < 4 && subgrid_geometry_edge_index(data.curlGeometry.data(), face, term) != edge)
				{
					++term;
				}
				if (term == 4)
				{
					throw std::logic_error("CartesianSubgrid: validated transpose lost its spatial edge slot");
				}
				high[k] = static_cast<float>(term);
				low[k] = 0;
			}
		}
		CartesianSubgridDeviceCurl packedTranspose(data.transposeStart, data.transposeColumn,
			high, low, data.magnetic.size());
		const std::size_t oldBytes = curl.DeviceBytes() + transpose.DeviceBytes();
		const std::size_t packedBytes = packedCurl.DeviceBytes() + packedTranspose.DeviceBytes() +
			data.curlGeometry.size() * sizeof(std::uint32_t) + data.curlMetrics.size() * sizeof(float);
		// A simple uniform grid may already have smaller translated patterns.
		if (packedBytes < oldBytes)
		{
			curl = std::move(packedCurl);
			transpose = std::move(packedTranspose);
			result.spatialRows = spatialRows;
			result.legacyBytes = oldBytes;
			result.geometry = data.curlGeometry;
			result.metrics = data.curlMetrics;
		}
	}

	result.curl = std::move(curl);
	result.transpose = std::move(transpose);
	return result;
}

//! Existing backend allocation/point operations; fields stay on that backend.
class CartesianSubgridExecution : public GPU_CartesianSubgrid
{
  public:
	CartesianSubgridExecution(GPU_Backend &backend, const CartesianSubgrid &grid)
		: CartesianSubgridExecution(backend, grid.DeviceStateData(), 0, 0, false, {}, &grid)
	{
	}

	CartesianSubgridExecution(GPU_Backend &backend, const CartesianSubgridDeviceData &data,
		std::size_t ownedElectric, std::size_t ownedMagnetic, bool domain,
		std::unique_ptr<GPU_CartesianSubgridRecurrence> recurrence = {},
		const CartesianSubgrid *streamingGrid = nullptr,
		const CartesianSubgridDeviceOperators *domainOperators = nullptr)
		: m_Backend(backend), m_Recurrence(std::move(recurrence)), m_Domain(domain),
		  m_OwnedElectric(domain ? ownedElectric : data.electric.size()),
		  m_OwnedMagnetic(domain ? ownedMagnetic : data.magnetic.size())
	{
		m_Backend.Activate();
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
		m_TrackSourceWork = m_PortCount != 0 || data.trackSourceWork;
		m_ElectricCount = data.electric.size();
		m_MagneticCount = data.magnetic.size();
		// Complete host field mirrors are allocated only when ReadFields is used.
		m_BlockCount = data.bandwidth.size();
		m_SourceCount = data.sourceEdge.size();
		const bool validStreamOwner = domain ? domainOperators && !streamingGrid : streamingGrid && !domainOperators;
		const bool stateOnly = data.curlColumn.empty() && data.curlWeight.empty() && data.curlLow.empty() &&
			data.curlStart == std::vector<std::uint32_t>{0} && data.transposeStart.empty() &&
			data.transposeColumn.empty() && data.transposeWeight.empty() && data.transposeLow.empty() &&
			data.curlGeometry.empty() && data.curlMetrics.empty();
		if ((domainOperators && !data.streamCurls) ||
			(data.streamCurls && (!validStreamOwner || !stateOnly)))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid state-only curl upload");
		}
		std::optional<CartesianSubgridDeviceOperators> prepared;
		if (!domainOperators)
		{
			prepared = data.streamCurls ? streamingGrid->DeviceCurls() : PackExpanded(data, domain);
		}
		const auto &operators = domainOperators ? *domainOperators : *prepared;
		const auto &curl = operators.curl;
		const auto &transpose = operators.transpose;
		m_SpatialRows = operators.spatialRows;
		if (m_SpatialRows)
		{
			m_Geometry = Upload(operators.geometry);
			m_Metrics = Upload(operators.metrics);
			const std::size_t bytes = curl.DeviceBytes() + transpose.DeviceBytes() +
				operators.geometry.size() * sizeof(std::uint32_t) + operators.metrics.size() * sizeof(float);
			std::cout << "CartesianSubgrid spatial device curl: rows=" << m_SpatialRows
				<< " metadata_bytes=" << bytes << " legacy_bytes=" << operators.legacyBytes << std::endl;
		}

		auto *electric = Upload(data.electric), *magnetic = Upload(data.magnetic);
		auto *increment = Upload(std::vector<float>(m_OwnedElectric, 0));
		m_Transpose = {Upload(transpose.start), Upload(transpose.offset), Upload(transpose.base),
			Upload(transpose.pattern), Upload(transpose.weight), Upload(transpose.residual),
			magnetic, increment, 1, 0, 0, m_Geometry, m_Metrics, 1};
		m_Curl = {Upload(curl.start), Upload(curl.offset), Upload(curl.base), Upload(curl.pattern),
			Upload(curl.weight), Upload(curl.residual), electric, magnetic, 0, 0, 1, m_Geometry, m_Metrics, 0};
		m_Mass = {Upload(data.blockStart), Upload(data.blockEdge), Upload(data.bandwidth),
			Upload(data.lowerStart), Upload(data.lower), increment};
		m_Update = {electric, increment, 0, 0};
		m_LedgerDevice = Upload(data.ledger);
		auto *sourceEdge = Upload(data.sourceEdge);
		auto *sourceCurl = Upload(std::vector<float>(data.sourceEdge.size(), 0));
		m_SourceCurl = {increment, sourceEdge, sourceCurl};
		m_Source = {electric, sourceEdge, Upload(data.sourceWeight), Upload(data.sourceOffStart),
			Upload(data.sourceOffEdge), Upload(data.sourceOffWeight), sourceCurl, m_LedgerDevice,
			0, 0, static_cast<std::uint32_t>(m_TrackSourceWork)};
		auto *portRHS = Upload(std::vector<float>(m_PortCount, 0));
		m_PortSolve = {Upload(data.portStart), Upload(data.portEdge), Upload(data.portWeight),
			Upload(data.rootConductance), Upload(data.portCurrent), Upload(data.portLower),
			electric, increment, portRHS, m_LedgerDevice, 0, 0, static_cast<std::uint32_t>(m_PortCount)};
		m_PortApply = {electric, Upload(data.responseStart), Upload(data.responseEdge),
			Upload(data.responseWeight), portRHS, 0};
		if (m_Recurrence)
		{
			auto owned = Fields();
			owned.electricCount = m_OwnedElectric;
			owned.magneticCount = m_OwnedMagnetic;
			m_Recurrence->Bind(m_Backend, owned, increment);
		}
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
		return {m_Update.electric, m_Curl.output, m_ElectricCount, m_MagneticCount,
			m_FixedStep != 0 ? m_FixedStep : m_PortStep};
	}

	void UpdateElectric(double timestep, double sourceAmplitude) override
	{
		if (m_Domain && m_OwnedMagnetic < m_MagneticCount && m_MagneticGhostEpoch != m_Epoch)
		{
			throw std::invalid_argument("CartesianSubgrid: stale magnetic ghost epoch");
		}
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
		m_Backend.Activate();
		m_Backend.RunPoints(GPU_OP_SUBGRID_CURL, &m_Transpose, 0, m_OwnedElectric, false);
		if (m_TrackSourceWork && m_SourceCount)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_SOURCE_CURL, &m_SourceCurl, 0, m_SourceCount, false);
		}
		if (m_Recurrence)
		{
			m_Recurrence->PrepareElectric(timestep);
		}
		m_Backend.RunPoints(GPU_OP_SUBGRID_MASS, &m_Mass, 0, m_BlockCount, false);
		if (m_PortCount)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_PORT_SOLVE, &m_PortSolve, 0, 1, false);
		}
		m_Backend.RunPoints(GPU_OP_SUBGRID_ELECTRIC, &m_Update, 0, m_OwnedElectric, false);
		if (m_PortCount)
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_PORT_APPLY, &m_PortApply, 0, m_PortCount, true);
		}
		if (m_Recurrence)
		{
			m_Recurrence->FinishElectric(timestep);
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
		if (m_Domain && m_OwnedElectric < m_ElectricCount && m_ElectricGhostEpoch != m_Epoch + 1)
		{
			throw std::invalid_argument("CartesianSubgrid: stale electric ghost epoch");
		}
		m_Backend.Activate();
		if (m_Recurrence)
		{
			m_Recurrence->UpdateMagnetic(m_Curl, timestep);
		}
		else
		{
			m_Backend.RunPoints(GPU_OP_SUBGRID_CURL, &m_Curl, 0, m_OwnedMagnetic, false);
		}
		m_PendingStep = 0;
		++m_Epoch;
	}

	void AcknowledgeGhosts(bool electric, std::size_t epoch) override
	{
		if (!m_Domain || (electric ? m_PendingStep == 0 || epoch != m_Epoch + 1 :
			m_PendingStep != 0 || epoch != m_Epoch))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid ghost completion phase/epoch");
		}
		(electric ? m_ElectricGhostEpoch : m_MagneticGhostEpoch) = epoch;
	}

	void ConfigureEnergy(const CartesianSubgrid &grid) override
	{
		if (m_Domain || grid.EdgeCount() != m_ElectricCount || grid.FaceCount() != m_MagneticCount)
		{
			throw std::invalid_argument("CartesianSubgrid: global energy grid differs from single-device fields");
		}
		std::vector<CartesianSubgridEnergyCoupling> terms;
		for (const auto &term : grid.ElectricEnergyCouplings())
		{
			if (term.first.second > std::numeric_limits<std::uint32_t>::max())
			{
				throw std::invalid_argument("CartesianSubgrid: energy pair exceeds indexing");
			}
			terms.push_back({static_cast<std::uint32_t>(term.first.first),
				static_cast<std::uint32_t>(term.first.second), term.second});
		}
		ConfigureEnergyCouplings(terms);
	}
	void ConfigureEnergyCouplings(const std::vector<CartesianSubgridEnergyCoupling> &terms) override
	{
		if (m_Started || m_PendingStep != 0 || m_EnergyAttempted)
		{
			throw std::logic_error("CartesianSubgrid: energy configuration is startup-only and single-use");
		}
		m_EnergyAttempted = true;
		auto owned = Fields();
		owned.electricCount = m_OwnedElectric;
		owned.magneticCount = m_OwnedMagnetic;
		m_Energy = GPU_CreateCartesianSubgridEnergy(m_Backend, owned, terms);
	}
	std::array<double, 2> ReadEnergy() override
	{
		if (!m_Energy || m_PendingStep != 0)
		{
			throw std::logic_error("CartesianSubgrid: field energy requires configured complete step");
		}
		if (m_Recurrence)
		{
			m_Recurrence->CheckObservation();
		}
		return m_Energy->Read();
	}
	std::uint64_t EnergyDownloadBytes() const override
	{
		return m_Energy ? m_Energy->DownloadBytes() : 0;
	}
	void ConfigureObservations(const std::vector<std::uint32_t> &electric,
		const std::vector<std::uint32_t> &magnetic) override
	{
		if (m_Started || m_PendingStep != 0 || m_ObservationAttempted)
		{
			throw std::logic_error("CartesianSubgrid: selective observation configuration requires unused startup");
		}
		const std::vector<std::uint32_t> *indices[2] = {&electric, &magnetic};
		const std::size_t limits[2] = {m_OwnedElectric, m_OwnedMagnetic};
		for (unsigned int field = 0; field < 2; ++field)
		{
			if (indices[field]->size() > std::numeric_limits<std::uint32_t>::max() ||
				std::any_of(indices[field]->begin(), indices[field]->end(),
					[&](std::uint32_t index) { return index >= limits[field]; }))
			{
				throw std::invalid_argument("CartesianSubgrid: selective observation index outside owned states");
			}
		}
		// Backend allocations live until backend destruction; a failed partial
		// configuration must not be retried as though that storage were freed.
		m_ObservationAttempted = true;
		m_Backend.Activate();
		for (unsigned int field = 0; field < 2; ++field)
		{
			m_ObservationCounts[field] = indices[field]->size();
			if (!indices[field]->empty())
			{
				m_ObservationGather[field] = {field ? m_Curl.output : m_Update.electric,
					Upload(*indices[field]), Upload(std::vector<float>(indices[field]->size(), 0))};
			}
		}
		m_ObservationReady = true;
	}
	void ReadObservations(std::vector<float> &electric, std::vector<float> &magnetic) override
	{
		if (!m_ObservationReady || m_PendingStep != 0 || &electric == &magnetic)
		{
			throw std::logic_error("CartesianSubgrid: selective observation requires configuration and complete step");
		}
		if (m_Recurrence)
		{
			m_Recurrence->CheckObservation();
		}
		std::vector<float> *output[2] = {&electric, &magnetic};
		m_Backend.Activate();
		for (unsigned int field = 0; field < 2; ++field)
		{
			const auto count = m_ObservationCounts[field];
			output[field]->resize(count);
			if (count)
			{
				auto &gather = m_ObservationGather[field];
				m_Backend.RunPoints(GPU_OP_SUBGRID_GATHER, &gather, 0, count, false);
				m_Backend.DeviceDownload(output[field]->data(), gather.curl, count * sizeof(float));
			}
		}
	}

	std::size_t AddObservationDFT(bool electric, const std::vector<std::uint32_t> &indices,
		std::size_t frequencies, std::size_t budget) override
	{
		if (m_Started || m_PendingStep != 0 || m_DFTStarted || m_DFTFailed)
		{
			throw std::logic_error("CartesianSubgrid: DFT registration requires unused startup");
		}
		const auto bytes = CartesianSubgridDFTBytes(indices.size(), frequencies);
		const auto owned = electric ? m_OwnedElectric : m_OwnedMagnetic;
		if (m_DFTBytes > budget || bytes > budget - m_DFTBytes ||
			std::any_of(indices.begin(), indices.end(), [&](std::uint32_t index) { return index >= owned; }))
		{
			throw std::invalid_argument("CartesianSubgrid: DFT budget/owned selection invalid");
		}
		try
		{
			auto accumulator = GPU_CreateCartesianSubgridDFT(m_Backend,
				electric ? m_Update.electric : m_Curl.output, owned, indices, frequencies);
			m_DFT.push_back({std::move(accumulator), frequencies, {}});
			m_DFTBytes += bytes;
		}
		catch (...)
		{
			// Allocations remain backend-owned. Do not retry a partial upload.
			m_DFTFailed = true;
			throw;
		}
		return m_DFT.size() - 1;
	}
	void AccumulateObservationDFT(std::size_t id, const std::vector<std::complex<float>> &weights) override
	{
		CheckDFT(id);
		auto &group = m_DFT[id];
		CartesianSubgridDFTWeights(weights, group.frequencies);
		if (group.lastEpoch && *group.lastEpoch == m_Epoch)
		{
			throw std::logic_error("CartesianSubgrid: duplicate DFT sample epoch");
		}
		try
		{
			group.accumulator->Accumulate(weights);
			group.lastEpoch = m_Epoch;
			m_DFTStarted = true;
		}
		catch (...)
		{
			m_DFTFailed = true;
			throw;
		}
	}
	void ReadObservationDFT(std::size_t id, std::size_t frequency,
		std::vector<std::complex<double>> &sums) override
	{
		CheckDFT(id);
		m_DFT[id].accumulator->Read(frequency, sums);
	}
	std::uint64_t ObservationDFTDownloadBytes() const override
	{
		std::uint64_t result = 0;
		for (const auto &group : m_DFT)
		{
			result += group.accumulator->DownloadBytes();
		}
		return result;
	}

	std::size_t HostObservationStorageBytes() const override
	{
		return (m_Electric.capacity() + m_Magnetic.capacity()) * sizeof(float);
	}
	void ReadFields(CartesianSubgrid &grid) override
	{
		if (m_Domain)
		{
			throw std::invalid_argument("CartesianSubgrid: domain observations require an owning coordinator");
		}
		if (m_PendingStep != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: device observation requires a complete timestep");
		}
		if (m_Recurrence)
		{
			m_Recurrence->CheckObservation();
		}
		m_Backend.Activate();
		m_Electric.resize(m_ElectricCount);
		m_Magnetic.resize(m_MagneticCount);
		m_Backend.DeviceDownload(m_Electric.data(), m_Update.electric, m_Electric.size() * sizeof(float));
		m_Backend.DeviceDownload(m_Magnetic.data(), m_Curl.output, m_Magnetic.size() * sizeof(float));
		grid.ImportDeviceFields(m_Electric.data(), m_Electric.size(), m_Magnetic.data(), m_Magnetic.size());
		if (m_TrackSourceWork)
		{
			const auto ledger = ReadLedger();
			grid.ImportDeviceLedger(ledger[0], ledger[1]);
		}
	}
	std::array<double, 2> ReadLedger() override
	{
		if (m_PendingStep != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: ledger requires a complete timestep");
		}
		m_Backend.Activate();
		m_Backend.DeviceDownload(m_Ledger.data(), m_LedgerDevice, sizeof(m_Ledger));
		return {m_Ledger[0], m_Ledger[2]};
	}

	std::size_t StructuredRows() const override
	{
		return m_SpatialRows;
	}

	std::size_t StorageBytes() const override
	{
		return m_StorageBytes + m_DFTBytes + (m_Recurrence ? m_Recurrence->StorageBytes() : 0) +
			(m_Energy ? m_Energy->StorageBytes() : 0);
	}
	std::array<double, 2> ReadAdvanceLedger() override
	{
		if (!m_Recurrence || m_PendingStep != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: recurrence diagnostics require a complete advanced step");
		}
		m_Recurrence->CheckObservation();
		m_Backend.Activate();
		return m_Recurrence->ReadLedger();
	}
	std::array<double, 2> ReadAdvanceSurfaceLedger() override
	{
		if (!m_Recurrence || m_PendingStep != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: surface diagnostics require a complete advanced step");
		}
		m_Recurrence->CheckObservation();
		m_Backend.Activate();
		return m_Recurrence->ReadSurfaceLedger();
	}

  private:

	void CheckDFT(std::size_t id) const
	{
		if (m_DFTFailed || id >= m_DFT.size() || m_PendingStep != 0)
		{
			throw std::logic_error("CartesianSubgrid: DFT requires a valid group and complete step");
		}
		if (m_Recurrence)
		{
			m_Recurrence->CheckObservation();
		}
	}
	struct DFTGroup
	{
		std::unique_ptr<GPU_CartesianSubgridDFT> accumulator;
		std::size_t frequencies;
		std::optional<std::size_t> lastEpoch;
	};
	std::vector<DFTGroup> m_DFT;
	std::size_t m_DFTBytes = 0;
	bool m_DFTStarted = false, m_DFTFailed = false;

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

	std::uint32_t *m_Geometry = nullptr;
	float *m_Metrics = nullptr;
	std::size_t m_SpatialRows = 0;
	GPU_Backend &m_Backend;
	std::unique_ptr<GPU_CartesianSubgridRecurrence> m_Recurrence;
	std::unique_ptr<GPU_CartesianSubgridEnergy> m_Energy;
	bool m_EnergyAttempted = false;
	std::size_t m_StorageBytes = 0, m_BlockCount = 0, m_SourceCount = 0, m_PortCount = 0;
	double m_StableStep = 0, m_PortStep = 0, m_FixedStep = 0, m_PendingStep = 0;
	bool m_Started = false, m_Domain = false, m_TrackSourceWork = false;
	std::size_t m_OwnedElectric = 0, m_OwnedMagnetic = 0, m_Epoch = 0;
	std::size_t m_ElectricCount = 0, m_MagneticCount = 0;
	std::size_t m_ElectricGhostEpoch = std::numeric_limits<std::size_t>::max();
	std::size_t m_MagneticGhostEpoch = std::numeric_limits<std::size_t>::max();
	std::vector<float> m_Electric, m_Magnetic;
	std::array<Op_SubgridGather, 2> m_ObservationGather{};
	std::array<std::size_t, 2> m_ObservationCounts{};
	bool m_ObservationAttempted = false, m_ObservationReady = false;
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

namespace
{
void ValidateDomain(const CartesianSubgridDomainData &domain, bool allowStream = false)
{
	const auto &data = domain.fields;
	const auto require = [](bool condition, const char *message)
	{
		if (!condition)
		{
			throw std::invalid_argument(message);
		}
	};
	require(domain.ownedElectric && domain.ownedElectric <= data.electric.size() &&
		domain.ownedMagnetic && domain.ownedMagnetic <= data.magnetic.size(),
		"CartesianSubgrid: invalid domain owned counts");
	require((data.streamCurls ? allowStream :
		(data.transposeStart.size() == domain.ownedElectric + 1 &&
			data.curlStart.size() == domain.ownedMagnetic + 1)) &&
		data.stableTimestep > 0 && std::isfinite(data.stableTimestep),
		"CartesianSubgrid: invalid domain rows/timestep");
	if (domain.operators)
	{
		const auto &operators = *domain.operators;
		require(data.streamCurls && operators.geometry.empty() && operators.metrics.empty() &&
			operators.spatialRows == 0, "CartesianSubgrid: invalid compressed domain upload mode");
		operators.curl.Validate(domain.ownedMagnetic, data.electric.size());
		operators.transpose.Validate(domain.ownedElectric, data.magnetic.size());
	}
	const std::size_t ports = data.rootConductance.size();
	require(ports <= 64 && data.portCurrent.size() == ports && data.portLower.size() == ports * ports &&
		(ports ? data.portTimestep > 0 && std::isfinite(data.portTimestep) &&
			data.portTimestep <= data.stableTimestep : data.portTimestep == 0),
		"CartesianSubgrid: invalid domain port counts/timestep");
	const auto sparse = [&](const std::vector<std::uint32_t> &start,
		const std::vector<std::uint32_t> &edge, const std::vector<float> &weight, std::size_t rows)
	{
		require(start.size() == rows + 1 && start.front() == 0 && start.back() == edge.size() &&
			edge.size() == weight.size() && std::is_sorted(start.begin(), start.end()),
			"CartesianSubgrid: malformed local port/source sparse rows");
		for (std::size_t index : edge)
		{
			require(index < domain.ownedElectric, "CartesianSubgrid: local port/source response writes a ghost");
		}
	};
	sparse(data.portStart, data.portEdge, data.portWeight, ports);
	sparse(data.responseStart, data.responseEdge, data.responseWeight, ports);
	for (std::size_t i = 0; i < ports; ++i)
	{
		require(data.rootConductance[i] > 0 && data.portLower[i * ports + i] > 0 &&
			data.portStart[i] < data.portStart[i + 1] && data.responseStart[i] < data.responseStart[i + 1],
			"CartesianSubgrid: nonpositive port factor or empty response");
		for (std::size_t j = i + 1; j < ports; ++j)
		{
			require(data.portLower[i * ports + j] == 0, "CartesianSubgrid: port upper triangle is not empty");
		}
	}
	require(data.ledger.size() == 4 && data.ledger[0] >= 0 &&
		(ports || data.trackSourceWork || data.ledger == std::vector<float>({0, 0, 0, 0})),
		"CartesianSubgrid: invalid/untracked domain ledger");
	require(data.blockStart.size() == data.bandwidth.size() + 1 &&
		data.lowerStart.size() == data.blockStart.size() && data.blockStart.front() == 0 &&
		data.lowerStart.front() == 0 && data.blockStart.back() == data.blockEdge.size() &&
		data.lowerStart.back() == data.lower.size(), "CartesianSubgrid: invalid domain factor bounds");
	std::vector<bool> seen(domain.ownedElectric, false);
	for (std::size_t block = 0; block < data.bandwidth.size(); ++block)
	{
		const std::size_t first = data.blockStart[block], stop = data.blockStart[block + 1];
		const std::size_t start = data.lowerStart[block], end = data.lowerStart[block + 1];
		require(first < stop && stop <= data.blockEdge.size() && start <= end && end <= data.lower.size(),
			"CartesianSubgrid: invalid domain factor range");
		const std::size_t entries = CartesianSubgridDetail::InterfaceFactorEntries(stop - first, data.bandwidth[block]);
		require(end - start == entries, "CartesianSubgrid: invalid domain factor size");
		for (std::size_t k = first; k < stop; ++k)
		{
			const std::size_t edge = data.blockEdge[k];
			require(edge < domain.ownedElectric && !seen[edge], "CartesianSubgrid: factor writes a ghost or overlapping component");
			seen[edge] = true;
			const float diagonal = data.lower[start + (k - first) * (data.bandwidth[block] + 1)];
			require(diagonal > 0, "CartesianSubgrid: nonpositive domain factor diagonal");
		}
	}
	require(data.sourceEdge.size() == data.sourceWeight.size(), "CartesianSubgrid: invalid domain sources");
	sparse(data.sourceOffStart, data.sourceOffEdge, data.sourceOffWeight, data.sourceEdge.size());
	for (std::size_t edge : data.sourceEdge)
	{
		require(edge < domain.ownedElectric, "CartesianSubgrid: source writes a ghost");
	}
	for (const auto *values : {&data.electric, &data.magnetic, &data.lower, &data.sourceWeight,
		&data.portWeight, &data.responseWeight, &data.rootConductance, &data.portCurrent,
		&data.portLower, &data.sourceOffWeight, &data.ledger})
	{
		for (float value : *values)
		{
			require(std::isfinite(value) && (value == 0 || std::abs(value) >= std::numeric_limits<float>::min()),
				"CartesianSubgrid: nonrepresentable domain state/factor/source");
		}
	}
}
}

std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridDomainExecution(
	GPU_Backend &backend, const CartesianSubgridDomainData &domain)
{
	ValidateDomain(domain, domain.operators.has_value());
	std::unique_ptr<GPU_CartesianSubgridRecurrence> recurrence;
	if (domain.loss)
	{
		if (!(domain.loss->timestep > 0) || domain.loss->timestep > domain.fields.stableTimestep)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid domain constitutive timestep");
		}
		recurrence = GPU_CreateCartesianSubgridConductivity(*domain.loss);
	}
	auto execution = std::make_unique<CartesianSubgridExecution>(backend, domain.fields,
		domain.ownedElectric, domain.ownedMagnetic, true, std::move(recurrence),
		nullptr, domain.operators ? &*domain.operators : nullptr);
	if (domain.loss)
	{
		execution->RequireFixedTimestep(domain.loss->timestep);
	}
	return execution;
}

std::unique_ptr<GPU_CartesianSubgrid> GPU_CreateCartesianSubgridAdvanced(
	GPU_Backend &backend, const CartesianSubgrid &grid, CartesianSubgridDeviceData data,
	std::unique_ptr<GPU_CartesianSubgridRecurrence> recurrence, double timestep)
{
	if (!recurrence || data.electric.size() != grid.EdgeCount() || data.magnetic.size() != grid.FaceCount() ||
		data.stableTimestep != grid.StableTimestep() || !(timestep > 0) || timestep > grid.StableTimestep())
	{
		throw std::invalid_argument("CartesianSubgrid: advanced execution geometry/recurrence mismatch");
	}
	CartesianSubgridDomainData checked;
	checked.ownedElectric = data.electric.size();
	checked.ownedMagnetic = data.magnetic.size();
	checked.fields = std::move(data);
	ValidateDomain(checked, true);
	auto execution = std::make_unique<CartesianSubgridExecution>(backend, checked.fields, 0, 0, false, std::move(recurrence), &grid);
	execution->RequireFixedTimestep(timestep);
	return execution;
}

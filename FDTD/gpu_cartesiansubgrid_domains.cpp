/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesian_subgrid_domains.h"
#include "cartesian_subgrid.h"
#include "gpu_backend.h"
#include "engine_plan.h"
#include "kernels/subgrid_update.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <limits>

namespace
{
using Field = CartesianSubgridField;
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}
std::size_t Slot(Field field)
{
	Require(field == Field::Electric || field == Field::Magnetic, "CartesianSubgrid: invalid exchange field");
	return field == Field::Electric ? 0 : 1;
}

const CartesianSubgridDomainData &ValidateWorkerDomain(GPU_Backend &backend, const CartesianSubgridDomainData &domain)
{
	for (Field field : {Field::Electric, Field::Magnetic})
	{
		const auto owned = field == Field::Electric ? domain.ownedElectric : domain.ownedMagnetic;
		const auto total = field == Field::Electric ? domain.fields.electric.size() : domain.fields.magnetic.size();
		Require(owned <= total, "CartesianSubgrid: invalid owned export count");
		const auto &indices = field == Field::Electric ? domain.sendElectric : domain.sendMagnetic;
		std::vector<bool> seen(owned, false);
		for (std::size_t index : indices)
		{
			Require(index < owned && !seen[index], "CartesianSubgrid: invalid owned export index");
			seen[index] = true;
		}
	}
	if (domain.pml)
	{
		Require(domain.pml->timestep > 0 && domain.pml->timestep <= domain.fields.stableTimestep &&
			(!domain.loss || domain.loss->timestep == domain.pml->timestep) &&
			(domain.fields.rootConductance.empty() || domain.fields.portTimestep == domain.pml->timestep),
			"CartesianSubgrid: invalid worker PML timestep");
		for (const auto &item : {std::make_pair(&domain.pml->electric, domain.ownedElectric),
			std::make_pair(&domain.pml->magnetic, domain.ownedMagnetic)})
		{
			Require(std::all_of(item.first->target.begin(), item.first->target.end(),
				[&](std::uint32_t target) { return target < item.second; }),
				"CartesianSubgrid: a ghost cannot own PML history");
		}
	}
	// Fail before allocating fields if this backend cannot update ghost storage.
	backend.Activate();
	backend.DeviceUpload(nullptr, nullptr, 0);
	return domain;
}

class LocalWorker : public CartesianSubgridDomainWorker
{
  public:
	LocalWorker(GPU_Backend &backend, const CartesianSubgridDomainData &domain)
		: m_Backend(backend), m_Execution(CreateCartesianSubgridDomainExecution(backend, ValidateWorkerDomain(backend, domain))),
		  m_TracksLedger(domain.fields.trackSourceWork || !domain.fields.rootConductance.empty())
	{
		m_RecurrenceStep = domain.loss ? domain.loss->timestep : 0;
		m_Owned = {domain.ownedElectric, domain.ownedMagnetic};
		if (domain.pml)
		{
			m_PMLStep = domain.pml->timestep;
			m_Execution->RequireFixedTimestep(m_PMLStep);
			m_PML = GPU_CreateCartesianSubgridPML(backend, m_Execution->Fields(), *domain.pml);
		}
		const auto fields = m_Execution->Fields();
		m_Total = {fields.electricCount, fields.magneticCount};
		for (Field field : {Field::Electric, Field::Magnetic})
		{
			const auto slot = Slot(field);
			const auto &indices = field == Field::Electric ? domain.sendElectric : domain.sendMagnetic;
			m_Exports[slot] = indices.size();
			m_Gather[slot] = {field == Field::Electric ? fields.electric : fields.magnetic,
				Upload(indices), Upload(std::vector<float>(indices.size(), 0))};
		}
	}

	void Electric(double timestep, double amplitude, std::size_t targetEpoch) override
	{
		Require(!m_Pending && targetEpoch == m_Epoch + 1, "CartesianSubgrid: worker electric epoch mismatch");
		Require(!m_PML || (timestep == m_PMLStep && std::isfinite(static_cast<float>(amplitude)) &&
			m_Epoch < std::numeric_limits<unsigned int>::max()), "CartesianSubgrid: invalid worker PML step");
		Hook(false, GPU_PRE);
		m_Execution->UpdateElectric(timestep, amplitude);
		Hook(false, GPU_POST);
		Hook(false, GPU_APPLY);
		m_Pending = timestep;
	}
	void Magnetic(double timestep, std::size_t targetEpoch) override
	{
		Require(m_Pending == timestep && m_Pending != 0 && targetEpoch == m_Epoch + 1,
			"CartesianSubgrid: worker magnetic epoch mismatch");
		Hook(true, GPU_PRE);
		m_Execution->UpdateMagnetic(timestep);
		Hook(true, GPU_POST);
		Hook(true, GPU_APPLY);
		m_Pending = 0;
		++m_Epoch;
	}
	void RequireFixedTimestep(double timestep) override { m_Execution->RequireFixedTimestep(timestep); }
	double RecurrenceTimestep() const override { return m_RecurrenceStep; }
	double PMLTimestep() const override { return m_PMLStep; }
	std::size_t Owned(Field field) const override { return m_Owned[Slot(field)]; }
	std::size_t Exports(Field field) const override { return m_Exports[Slot(field)]; }
	std::size_t Ghosts(Field field) const override { return m_Total[Slot(field)] - Owned(field); }

	void ReadExports(Field field, std::size_t epoch, float *values, std::size_t count) override
	{
		CheckExchange(field, epoch);
		Require(count == Exports(field) && (!count || values), "CartesianSubgrid: invalid export buffer");
		if (count)
		{
			m_Backend.Activate();
			auto &gather = m_Gather[Slot(field)];
			m_Backend.RunPoints(GPU_OP_SUBGRID_GATHER, &gather, 0, count, false);
			m_Backend.DeviceDownload(values, gather.curl, count * sizeof(float));
		}
	}
	void InstallGhosts(Field field, std::size_t epoch, const float *values, std::size_t count) override
	{
		CheckExchange(field, epoch);
		Require(count == Ghosts(field) && (!count || values), "CartesianSubgrid: invalid ghost buffer");
		if (count)
		{
			Require(std::all_of(values, values + count, [](float value) { return std::isfinite(value); }),
				"CartesianSubgrid: nonfinite ghost packet");
			const auto fields = m_Execution->Fields();
			auto *base = field == Field::Electric ? fields.electric : fields.magnetic;
			m_Backend.Activate();
			m_Backend.DeviceUpload(base + Owned(field), values, count * sizeof(float));
		}
		m_Execution->AcknowledgeGhosts(field == Field::Electric, epoch);
	}
	void ReadOwned(std::size_t epoch, float *electric, std::size_t electricCount,
		float *magnetic, std::size_t magneticCount) override
	{
		Require(!m_Pending && epoch == m_Epoch, "CartesianSubgrid: mixed worker observation epoch");
		Require(electricCount == Owned(Field::Electric) && magneticCount == Owned(Field::Magnetic) &&
			electric && magnetic, "CartesianSubgrid: invalid owned observation buffers");
		const auto fields = m_Execution->Fields();
		m_Backend.Activate();
		m_Backend.DeviceDownload(electric, fields.electric, electricCount * sizeof(float));
		m_Backend.DeviceDownload(magnetic, fields.magnetic, magneticCount * sizeof(float));
	}
	std::size_t StorageBytes() const override
	{
		return m_StorageBytes + m_Execution->StorageBytes() + (m_PML ? m_PML->StorageBytes() : 0);
	}
	double ReadBoundaryWorkOut(std::size_t epoch) override
	{
		Require(m_PML && !m_Pending && epoch == m_Epoch, "CartesianSubgrid: mixed worker PML epoch");
		m_Backend.Activate();
		return m_PML->ReadBoundaryWorkOut();
	}
	std::array<double, 2> ReadLedger(std::size_t epoch) override
	{
		Require(!m_Pending && epoch == m_Epoch, "CartesianSubgrid: mixed worker ledger epoch");
		return m_TracksLedger ? m_Execution->ReadLedger() : std::array<double, 2>{0, 0};
	}
	std::array<double, 2> ReadAdvanceLedger(std::size_t epoch) override
	{
		Require(!m_Pending && epoch == m_Epoch, "CartesianSubgrid: mixed worker recurrence epoch");
		return m_Execution->ReadAdvanceLedger();
	}
	std::array<double, 2> ReadAdvanceSurfaceLedger(std::size_t epoch) override
	{
		Require(!m_Pending && epoch == m_Epoch, "CartesianSubgrid: mixed worker surface epoch");
		return m_Execution->ReadAdvanceSurfaceLedger();
	}

  private:
	void Hook(bool magnetic, GPU_Hook phase)
	{
		if (m_PML)
		{
			m_Backend.Activate();
			m_PML->Hook(magnetic, phase, Engine_Box{}, static_cast<unsigned int>(m_Epoch));
		}
	}
	void CheckExchange(Field field, std::size_t epoch) const
	{
		Slot(field);
		Require(field == Field::Electric ? m_Pending != 0 && epoch == m_Epoch + 1 :
			m_Pending == 0 && epoch == m_Epoch, "CartesianSubgrid: worker exchange phase/epoch mismatch");
	}
	template <class T> T *Upload(const std::vector<T> &values)
	{
		const std::size_t bytes = std::max(values.size() * sizeof(T), std::size_t(4));
		auto *buffer = m_Backend.DeviceAlloc(bytes, values.empty() ? nullptr : values.data());
		if (!buffer)
		{
			throw std::runtime_error("CartesianSubgrid: domain exchange allocation failed");
		}
		m_StorageBytes += bytes;
		return static_cast<T*>(buffer);
	}
	GPU_Backend &m_Backend;
	std::unique_ptr<GPU_CartesianSubgrid> m_Execution;
	std::unique_ptr<GPU_CartesianSubgridPML> m_PML;
	std::array<std::size_t, 2> m_Owned{}, m_Total{}, m_Exports{};
	std::array<Op_SubgridGather, 2> m_Gather{};
	std::size_t m_Epoch = 0, m_StorageBytes = 0;
	double m_Pending = 0, m_RecurrenceStep = 0, m_PMLStep = 0;
	bool m_TracksLedger = false;
};
}

std::unique_ptr<CartesianSubgridDomainWorker> CreateCartesianSubgridWorker(
	GPU_Backend &backend, const CartesianSubgridDomainData &domain)
{
	return std::make_unique<LocalWorker>(backend, domain);
}
std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridDomains(
	std::shared_ptr<const CartesianSubgridDomainPlan> plan, const std::vector<GPU_Backend*> &backends)
{
	Require(bool(plan) && backends.size() == plan->Domains().size(), "CartesianSubgrid: invalid domain backends");
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers;
	for (std::size_t k = 0; k < backends.size(); ++k)
	{
		Require(backends[k] != nullptr, "CartesianSubgrid: null domain backend");
		workers.push_back(CreateCartesianSubgridWorker(*backends[k], plan->Domains()[k]));
	}
	auto transport = CreateCartesianSubgridHostTransport(plan, workers);
	return CreateCartesianSubgridCoordinator(std::move(plan), std::move(workers), std::move(transport));
}

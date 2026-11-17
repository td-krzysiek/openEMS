/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesian_subgrid_domains.h"
#include "cartesian_subgrid.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

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

void CheckWorkers(const CartesianSubgridDomainPlan &plan,
	const std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> &workers)
{
	Require(workers.size() == plan.Domains().size(), "CartesianSubgrid: worker count differs from ownership plan");
	for (std::size_t k = 0; k < workers.size(); ++k)
	{
		Require(bool(workers[k]), "CartesianSubgrid: missing domain worker");
		const auto &d = plan.Domains()[k];
		Require(workers[k]->RecurrenceTimestep() == (d.loss ? d.loss->timestep : 0),
			"CartesianSubgrid: worker constitutive preparation differs from ownership plan");
		Require(workers[k]->PMLTimestep() == (d.pml ? d.pml->timestep : 0),
			"CartesianSubgrid: worker PML preparation differs from ownership plan");
		Require(workers[k]->Owned(Field::Electric) == d.ownedElectric &&
			workers[k]->Owned(Field::Magnetic) == d.ownedMagnetic &&
			workers[k]->Ghosts(Field::Electric) == d.globalElectric.size() - d.ownedElectric &&
			workers[k]->Ghosts(Field::Magnetic) == d.globalMagnetic.size() - d.ownedMagnetic &&
			workers[k]->Exports(Field::Electric) == d.sendElectric.size() &&
			workers[k]->Exports(Field::Magnetic) == d.sendMagnetic.size(),
			"CartesianSubgrid: worker storage differs from ownership plan");
	}
}

class HostTransport : public CartesianSubgridDomainTransport
{
  public:
	HostTransport(std::shared_ptr<const CartesianSubgridDomainPlan> plan,
		const std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> &workers)
		: m_Plan(std::move(plan))
	{
		Require(bool(m_Plan), "CartesianSubgrid: missing exchange plan");
		CheckWorkers(*m_Plan, workers);
		for (const auto &worker : workers)
		{
			m_Workers.push_back(worker.get());
		}
		for (Field field : {Field::Electric, Field::Magnetic})
		{
			auto &buffers = m_Buffers[Slot(field)];
			buffers.send.resize(workers.size());
			buffers.receive.resize(workers.size());
			for (std::size_t k = 0; k < workers.size(); ++k)
			{
				buffers.send[k].resize(workers[k]->Exports(field));
				buffers.receive[k].resize(workers[k]->Ghosts(field));
				buffers.bytes += (buffers.send[k].size() + buffers.receive[k].size()) * sizeof(float);
			}
		}
	}
	void Exchange(Field field, std::size_t epoch) override
	{
		auto &buffers = m_Buffers[Slot(field)];
		for (std::size_t k = 0; k < m_Workers.size(); ++k)
		{
			m_Workers[k]->ReadExports(field, epoch, buffers.send[k].data(), buffers.send[k].size());
		}
		for (const auto &transfer : m_Plan->Transfers(field))
		{
			buffers.receive[transfer.receiver][transfer.receiveOffset] = buffers.send[transfer.sender][transfer.sendOffset];
		}
		for (std::size_t k = 0; k < m_Workers.size(); ++k)
		{
			m_Workers[k]->InstallGhosts(field, epoch, buffers.receive[k].data(), buffers.receive[k].size());
		}
		m_Traffic += buffers.bytes;
	}
	std::uint64_t HostTrafficBytes() const override { return m_Traffic; }

  private:
	struct Buffers
	{
		std::vector<std::vector<float>> send, receive;
		std::uint64_t bytes = 0;
	};
	std::shared_ptr<const CartesianSubgridDomainPlan> m_Plan;
	std::vector<CartesianSubgridDomainWorker*> m_Workers;
	std::array<Buffers, 2> m_Buffers;
	std::uint64_t m_Traffic = 0;
};

//! Scheduling depends only on workers/transport, not local GPU handles.
class Coordinator : public GPU_CartesianSubgrid
{
  public:
	Coordinator(std::shared_ptr<const CartesianSubgridDomainPlan> plan,
		std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers,
		std::unique_ptr<CartesianSubgridDomainTransport> transport)
		: m_Plan(std::move(plan)), m_Workers(std::move(workers)), m_Transport(std::move(transport))
	{
		Require(bool(m_Plan) && bool(m_Transport), "CartesianSubgrid: missing coordinator plan/transport");
		CheckWorkers(*m_Plan, m_Workers);
		std::size_t maxE = 0, maxH = 0;
		for (const auto &worker : m_Workers)
		{
			maxE = std::max(maxE, worker->Owned(Field::Electric));
			maxH = std::max(maxH, worker->Owned(Field::Magnetic));
		}
		m_ReadElectric.resize(maxE);
		m_ReadMagnetic.resize(maxH);
		m_Electric.resize(m_Plan->ElectricCount());
		m_Magnetic.resize(m_Plan->MagneticCount());
		if (m_Plan->FactorTails().terms || m_Plan->StateTails().terms)
		{
			std::cout << "CartesianSubgrid domain plan FP32 FTZ tails: factors=" << m_Plan->FactorTails().terms
				<< ", states=" << m_Plan->StateTails().terms << std::endl;
		}
	}
	void Step(double timestep, double amplitude) override
	{
		UpdateElectric(timestep, amplitude);
		UpdateMagnetic(timestep);
	}
	void RequireFixedTimestep(double timestep) override
	{
		Healthy();
		CheckStep(timestep, 0);
		Require(m_Epoch == 0 && !m_Pending, "CartesianSubgrid: fixed coordinator timestep requires startup");
		try
		{
			for (auto &worker : m_Workers)
			{
				worker->RequireFixedTimestep(timestep);
			}
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
		m_Fixed = timestep;
	}
	void UpdateElectric(double timestep, double amplitude) override
	{
		Healthy();
		Require(!m_Pending, "CartesianSubgrid: coordinator requires the previous complete step");
		CheckStep(timestep, amplitude);
		try
		{
			m_Transport->Exchange(Field::Magnetic, m_Epoch);
			for (auto &worker : m_Workers)
			{
				worker->Electric(timestep, amplitude, m_Epoch + 1);
			}
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
		m_Pending = timestep;
	}
	void UpdateMagnetic(double timestep) override
	{
		Healthy();
		Require(m_Pending != 0 && m_Pending == timestep, "CartesianSubgrid: coordinator magnetic half-step mismatch");
		try
		{
			m_Transport->Exchange(Field::Electric, m_Epoch + 1);
			for (auto &worker : m_Workers)
			{
				worker->Magnetic(timestep, m_Epoch + 1);
			}
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
		m_Pending = 0;
		++m_Epoch;
	}
	void ReadFields(CartesianSubgrid &grid) override
	{
		Healthy();
		Require(!m_Pending, "CartesianSubgrid: coordinator observation requires a complete step");
		Require(grid.EdgeCount() == m_Electric.size() && grid.FaceCount() == m_Magnetic.size(),
			"CartesianSubgrid: coordinator observation dimensions differ");
		try
		{
			for (std::size_t k = 0; k < m_Workers.size(); ++k)
			{
				const auto &domain = m_Plan->Domains()[k];
				m_Workers[k]->ReadOwned(m_Epoch, m_ReadElectric.data(), domain.ownedElectric,
					m_ReadMagnetic.data(), domain.ownedMagnetic);
				for (std::size_t i = 0; i < domain.ownedElectric; ++i)
				{
					m_Electric[domain.globalElectric[i]] = m_ReadElectric[i];
				}
				for (std::size_t i = 0; i < domain.ownedMagnetic; ++i)
				{
					m_Magnetic[domain.globalMagnetic[i]] = m_ReadMagnetic[i];
				}
			}
			const auto ledger = ReadLedger();
			grid.ImportDeviceFields(m_Electric.data(), m_Electric.size(), m_Magnetic.data(), m_Magnetic.size());
			grid.ImportDeviceLedger(ledger[0], ledger[1]);
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
	}
	std::array<double, 2> ReadLedger() override
	{
		Healthy();
		Require(!m_Pending, "CartesianSubgrid: coordinator ledger requires a complete step");
		std::array<double, 2> result{};
		try
		{
			for (auto &worker : m_Workers)
			{
				const auto ledger = worker->ReadLedger(m_Epoch);
				Require(std::isfinite(ledger[0]) && ledger[0] >= 0 && std::isfinite(ledger[1]),
					"CartesianSubgrid: invalid worker work/loss ledger");
				result[0] += ledger[0];
				result[1] += ledger[1];
			}
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
		return result;
	}
	std::size_t StorageBytes() const override
	{
		std::size_t bytes = 0;
		for (const auto &worker : m_Workers)
		{
			bytes += worker->StorageBytes();
		}
		return bytes;
	}
	std::array<double, 2> ReadAdvanceLedger() override
	{
		return ReadConstitutiveLedger(false);
	}
	std::array<double, 2> ReadAdvanceSurfaceLedger() override
	{
		return ReadConstitutiveLedger(true);
	}
	double ReadBoundaryWorkOut() override
	{
		Healthy();
		Require(!m_Pending && m_Plan->PMLTimestep() > 0,
			"CartesianSubgrid: boundary observation requires a complete PML step");
		double result = 0;
		try
		{
			for (auto &worker : m_Workers)
			{
				const double value = worker->ReadBoundaryWorkOut(m_Epoch);
				Require(std::isfinite(value), "CartesianSubgrid: invalid worker PML work");
				result += value;
				Require(std::isfinite(result), "CartesianSubgrid: PML work overflow");
			}
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
		return result;
	}
	GPU_CartesianSubgridFields Fields() const override
	{
		throw std::logic_error("CartesianSubgrid: domains have no single borrowed field buffer; use the coordinator");
	}

  private:
	std::array<double, 2> ReadConstitutiveLedger(bool surface)
	{
		Healthy();
		Require(!m_Pending && m_Plan->RecurrenceTimestep() > 0,
			"CartesianSubgrid: recurrence observation requires a complete lossy step");
		std::array<double, 2> result{};
		try
		{
			for (auto &worker : m_Workers)
			{
				const auto value = surface ? worker->ReadAdvanceSurfaceLedger(m_Epoch) :
					worker->ReadAdvanceLedger(m_Epoch);
				for (std::size_t k = 0; k < result.size(); ++k)
				{
					Require(value[k] >= 0 && std::isfinite(value[k]),
						"CartesianSubgrid: invalid worker constitutive ledger");
					result[k] += value[k];
					Require(std::isfinite(result[k]), "CartesianSubgrid: constitutive ledger overflow");
				}
			}
		}
		catch (...)
		{
			m_Failed = true;
			throw;
		}
		return result;
	}
	void Healthy() const { Require(!m_Failed, "CartesianSubgrid: failed coordinator cannot continue"); }
	void CheckStep(double timestep, double amplitude) const
	{
		const float dt = static_cast<float>(timestep), low = static_cast<float>(timestep - double(dt));
		Require(timestep > 0 && timestep <= m_Plan->StableTimestep() && std::isfinite(dt) &&
			dt >= std::numeric_limits<float>::min() && std::isfinite(static_cast<float>(amplitude)) &&
			(timestep == double(dt) || std::abs(low) >= std::numeric_limits<float>::min()) &&
			(m_Plan->PortTimestep() == 0 || timestep == m_Plan->PortTimestep()) &&
			(m_Plan->RecurrenceTimestep() == 0 || timestep == m_Plan->RecurrenceTimestep()) &&
			(m_Plan->PMLTimestep() == 0 || timestep == m_Plan->PMLTimestep()) &&
			(m_Fixed == 0 || m_Fixed == timestep), "CartesianSubgrid: invalid coordinator timestep/source");
	}
	std::shared_ptr<const CartesianSubgridDomainPlan> m_Plan;
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> m_Workers;
	std::unique_ptr<CartesianSubgridDomainTransport> m_Transport;
	std::vector<float> m_Electric, m_Magnetic, m_ReadElectric, m_ReadMagnetic;
	std::size_t m_Epoch = 0;
	double m_Pending = 0, m_Fixed = 0;
	bool m_Failed = false;
};
}

std::unique_ptr<CartesianSubgridDomainTransport> CreateCartesianSubgridHostTransport(
	std::shared_ptr<const CartesianSubgridDomainPlan> plan,
	const std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> &workers)
{
	return std::make_unique<HostTransport>(std::move(plan), workers);
}
std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridCoordinator(
	std::shared_ptr<const CartesianSubgridDomainPlan> plan,
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers,
	std::unique_ptr<CartesianSubgridDomainTransport> transport)
{
	return std::make_unique<Coordinator>(std::move(plan), std::move(workers), std::move(transport));
}

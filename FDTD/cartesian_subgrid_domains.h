/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_DOMAINS_H
#define CARTESIAN_SUBGRID_DOMAINS_H

#include "gpu_cartesiansubgrid.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <memory>
#include <vector>

//! Compact ownership is independent of backend, process and host placement.
enum class CartesianSubgridField { Electric, Magnetic };

struct CartesianSubgridDomainData
{
	CartesianSubgridDeviceData fields;
	std::size_t ownedElectric = 0, ownedMagnetic = 0;
	//! Owned entries first, then read-only ghosts, in each local field array.
	std::vector<std::uint32_t> globalElectric, globalMagnetic;
	//! Unique owned local indices exported to other workers.
	std::vector<std::uint32_t> sendElectric, sendMagnetic;
};

//! Immutable fixed ownership/dependency plan; not a load-balancing optimizer.
class CartesianSubgridDomainPlan
{
  public:
	struct Transfer
	{
		std::uint32_t sender, sendOffset, receiver, receiveOffset;
	};
	//! Cuts are assignment planes, not new mesh faces. Leaves stay indivisible.
	CartesianSubgridDomainPlan(const CartesianSubgrid &grid, unsigned int axis,
		const std::vector<double> &cuts);
	const std::vector<CartesianSubgridDomainData> &Domains() const { return m_Domains; }
	const std::vector<Transfer> &Transfers(CartesianSubgridField field) const
	{
		if (field != CartesianSubgridField::Electric && field != CartesianSubgridField::Magnetic)
		{
			throw std::invalid_argument("CartesianSubgrid: invalid dependency field");
		}
		return field == CartesianSubgridField::Electric ? m_ElectricTransfers : m_MagneticTransfers;
	}
	std::size_t ElectricCount() const { return m_ElectricCount; }
	std::size_t MagneticCount() const { return m_MagneticCount; }
	std::size_t ReassignedElectricCount() const { return m_ReassignedElectric; }
	double StableTimestep() const { return m_StableStep; }
	double PortTimestep() const { return m_PortStep; }
	const CartesianSubgridFloat::TailInfo &FactorTails() const { return m_FactorTails; }
	const CartesianSubgridFloat::TailInfo &StateTails() const { return m_StateTails; }

  private:
	std::vector<CartesianSubgridDomainData> m_Domains;
	std::vector<Transfer> m_ElectricTransfers, m_MagneticTransfers;
	std::size_t m_ElectricCount = 0, m_MagneticCount = 0, m_ReassignedElectric = 0;
	double m_StableStep = 0, m_PortStep = 0;
	CartesianSubgridFloat::TailInfo m_FactorTails, m_StateTails;
};

//! Ordered half-step submission; a remote adapter may implement this contract.
//! Export/owned reads and ghost installation complete before returning. A queued
//! epoch is not a physical-device fence; processing must wait for owned reads.
//! Ports keep their complete coupled solve with one owner; other histories reject.
class CartesianSubgridDomainWorker
{
  public:
	virtual ~CartesianSubgridDomainWorker() = default;
	virtual void Electric(double timestep, double amplitude, std::size_t targetEpoch) = 0;
	virtual void Magnetic(double timestep, std::size_t targetEpoch) = 0;
	virtual void RequireFixedTimestep(double timestep) = 0;
	virtual std::size_t Owned(CartesianSubgridField field) const = 0;
	virtual std::size_t Exports(CartesianSubgridField field) const = 0;
	virtual std::size_t Ghosts(CartesianSubgridField field) const = 0;
	virtual void ReadExports(CartesianSubgridField field, std::size_t epoch, float *values, std::size_t count) = 0;
	virtual void InstallGhosts(CartesianSubgridField field, std::size_t epoch, const float *values, std::size_t count) = 0;
	virtual void ReadOwned(std::size_t epoch, float *electric, std::size_t electricCount,
		float *magnetic, std::size_t magneticCount) = 0;
	virtual std::array<double, 2> ReadLedger(std::size_t epoch) = 0;
	virtual std::size_t StorageBytes() const = 0;
};

//! Scheduling and ownership do not depend on the transport implementation.
class CartesianSubgridDomainTransport
{
  public:
	virtual ~CartesianSubgridDomainTransport() = default;
	virtual void Exchange(CartesianSubgridField field, std::size_t epoch) = 0;
	virtual std::uint64_t HostTrafficBytes() const = 0;
};

//! Backends outlive their workers. A repeated backend pointer means logical domains.
std::unique_ptr<CartesianSubgridDomainWorker> CreateCartesianSubgridWorker(
	GPU_Backend &backend, const CartesianSubgridDomainData &domain);
std::unique_ptr<CartesianSubgridDomainTransport> CreateCartesianSubgridHostTransport(
	std::shared_ptr<const CartesianSubgridDomainPlan> plan,
	const std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> &workers);
std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridCoordinator(
	std::shared_ptr<const CartesianSubgridDomainPlan> plan,
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers,
	std::unique_ptr<CartesianSubgridDomainTransport> transport);
//! Local convenience factory; no native XML/API or distributed extension selection yet.
std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridDomains(
	std::shared_ptr<const CartesianSubgridDomainPlan> plan, const std::vector<GPU_Backend*> &backends);

//! Internal local-domain execution; validates owned rows/factors before upload.
std::unique_ptr<GPU_CartesianSubgrid> CreateCartesianSubgridDomainExecution(
	GPU_Backend &backend, const CartesianSubgridDomainData &domain);
#endif

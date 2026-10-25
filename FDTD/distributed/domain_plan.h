// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_DOMAIN_PLAN_H
#define OPENEMS_DISTRIBUTED_DOMAIN_PLAN_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace openems
{
namespace distributed
{

size_t CheckedAdd(size_t a, size_t b);
size_t CheckedMultiply(size_t a, size_t b);

enum class FieldPhase : uint64_t
{
	Voltage = 0,
	Current = 1
};

struct GridShape
{
	uint32_t nx, ny, nz;
	size_t Plane() const;
	size_t Nodes() const;
	void Validate() const;
};

// X ranges own native grid nodes, not physical cell intervals. Current updates
// additionally exclude the true global last line in every direction.
struct Domain
{
	uint32_t id, begin, end;
	GridShape grid;
	size_t OwnedNodes() const;
	size_t StorageNodes() const;
	size_t LocalIndex(uint32_t component, uint32_t local_x, uint32_t y, uint32_t z) const;
	uint32_t ToLocal(uint32_t global_x) const;
	bool Owns(uint32_t global_x) const;
};

struct HaloRecipe
{
	uint32_t channel, source, destination;
	FieldPhase phase;
	uint32_t source_x, destination_x; // local storage indices including ghosts
	size_t plane_samples;
	size_t Samples() const;
};

class PartitionPlan
{
public:
	PartitionPlan(uint64_t id, GridShape shape, const std::vector<uint32_t>& cuts);

	uint64_t Id() const
	{
		return m_Id;
	}

	const GridShape& Shape() const
	{
		return m_Shape;
	}

	const std::vector<Domain>& Domains() const
	{
		return m_Domains;
	}

	const std::vector<HaloRecipe>& Halos() const
	{
		return m_Halos;
	}

	uint32_t Owner(uint32_t global_x) const;

private:
	uint64_t m_Id;
	GridShape m_Shape;
	std::vector<Domain> m_Domains;
	std::vector<HaloRecipe> m_Halos;
};

struct Placement
{
	uint32_t domain, process;
	std::string host, backend, resource;
	size_t memory_budget_bytes, auxiliary_bytes;
};

// Conservative uncompressed baseline, not a prediction for opaque/PML/DFT state.
size_t BasicDomainBytes(const Domain& domain);
void ValidatePlacement(const PartitionPlan& plan, const std::vector<Placement>& placement);

struct ExchangeEpoch
{
	uint64_t plan, step;
	FieldPhase phase;
	bool operator==(const ExchangeEpoch& other) const;
};

struct MessageId
{
	ExchangeEpoch epoch;
	uint32_t channel, source, destination;
	bool operator<(const MessageId& other) const;
	bool operator==(const MessageId& other) const;
};

}
}
#endif

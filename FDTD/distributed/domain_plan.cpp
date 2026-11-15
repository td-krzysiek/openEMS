// SPDX-License-Identifier: GPL-3.0-or-later
#include "domain_plan.h"
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

namespace openems
{
namespace distributed
{

size_t CheckedAdd(size_t a, size_t b)
{
	if (b > std::numeric_limits<size_t>::max() - a)
	{
		throw std::overflow_error("distributed size addition overflow");
	}
	return a + b;
}

size_t CheckedMultiply(size_t a, size_t b)
{
	if (a && b > std::numeric_limits<size_t>::max() / a)
	{
		throw std::overflow_error("distributed size multiplication overflow");
	}
	return a * b;
}

size_t GridShape::Plane() const
{
	return CheckedMultiply(ny, nz);
}

size_t GridShape::Nodes() const
{
	return CheckedMultiply(nx, Plane());
}

void GridShape::Validate() const
{
	if (nx < 2 || ny < 2 || nz < 2)
	{
		throw std::invalid_argument("Cartesian domains require at least two grid lines per axis");
	}
	CheckedMultiply(Nodes(), 18 * sizeof(float));
}

size_t Domain::OwnedNodes() const
{
	return CheckedMultiply(end - begin, grid.Plane());
}

size_t Domain::StorageNodes() const
{
	return CheckedMultiply(CheckedAdd(end - begin, 2), grid.Plane());
}

bool Domain::Owns(uint32_t x) const
{
	return x >= begin && x < end;
}

uint32_t Domain::ToLocal(uint32_t x) const
{
	if (!Owns(x))
	{
		throw std::out_of_range("sample is not owned by this domain");
	}
	return x - begin + 1;
}

size_t Domain::LocalIndex(uint32_t c, uint32_t x, uint32_t y, uint32_t z) const
{
	if (c > 2 || size_t(x) > size_t(end - begin) + 1 || y >= grid.ny || z >= grid.nz)
	{
		throw std::out_of_range("local field index outside domain storage");
	}
	return CheckedAdd(CheckedMultiply(c, StorageNodes()),
	    CheckedAdd(CheckedMultiply(x, grid.Plane()), CheckedAdd(CheckedMultiply(y, grid.nz), z)));
}

size_t HaloRecipe::Samples() const
{
	return CheckedMultiply(2, plane_samples);
}

PartitionPlan::PartitionPlan(uint64_t id, GridShape shape, const std::vector<uint32_t>& cuts)
    : m_Id(id), m_Shape(shape)
{
	shape.Validate();
	if (!id || cuts.size() < 2 || cuts.front() != 0 || cuts.back() != shape.nx)
	{
		throw std::invalid_argument("partition must have an ID and cover all global X lines");
	}
	if (cuts.size() - 1 > std::numeric_limits<uint32_t>::max() / 2)
	{
		throw std::overflow_error("too many distributed domains");
	}
	for (size_t i = 0; i + 1 < cuts.size(); ++i)
	{
		if (cuts[i] >= cuts[i + 1] || cuts[i + 1] > shape.nx ||
		    uint64_t(cuts[i + 1]) - cuts[i] + 2 > UINT32_MAX)
		{
			throw std::invalid_argument(
			    "partition cuts must be strictly ordered with nonempty storage");
		}
		Domain domain = { uint32_t(i), cuts[i], cuts[i + 1], shape };
		BasicDomainBytes(domain);
		m_Domains.push_back(domain);
	}
	for (uint32_t i = 0; i + 1 < m_Domains.size(); ++i)
	{
		// Current HIP curl: V reads backward I; I reads forward V in X.
		HaloRecipe v = { 2 * i,
			i + 1,
			i,
			FieldPhase::Voltage,
			1,
			m_Domains[i].end - m_Domains[i].begin + 1,
			shape.Plane() };
		HaloRecipe h = { 2 * i + 1,
			i,
			i + 1,
			FieldPhase::Current,
			m_Domains[i].end - m_Domains[i].begin,
			0,
			shape.Plane() };
		m_Halos.push_back(v);
		m_Halos.push_back(h);
	}
}

uint32_t PartitionPlan::Owner(uint32_t x) const
{
	if (x >= m_Shape.nx)
	{
		throw std::out_of_range("global X line outside plan");
	}
	const auto it = std::upper_bound(m_Domains.begin(),
	    m_Domains.end(),
	    x,
	    [](uint32_t value, const Domain& d)
	{
		return value < d.end;
	});
	return it->id;
}

size_t BasicDomainBytes(const Domain& d)
{
	return CheckedAdd(CheckedMultiply(d.StorageNodes(), 6 * sizeof(float)),
	    CheckedMultiply(d.OwnedNodes(), 12 * sizeof(float)));
}

void ValidatePlacement(const PartitionPlan& plan, const std::vector<Placement>& placement)
{
	if (placement.size() != plan.Domains().size())
	{
		throw std::invalid_argument("placement must name every domain exactly once");
	}
	std::set<uint32_t> seen;
	std::map<std::pair<std::string, std::string>, std::pair<size_t, size_t>> resources;
	for (const auto& p : placement)
	{
		if (p.domain >= plan.Domains().size() || !seen.insert(p.domain).second || p.host.empty() ||
		    p.backend.empty() || p.resource.empty() || !p.memory_budget_bytes)
		{
			throw std::invalid_argument("invalid or duplicate domain placement");
		}
		auto& resource = resources[std::make_pair(p.host, p.resource)];
		if (resource.second && resource.second != p.memory_budget_bytes)
		{
			throw std::invalid_argument("inconsistent resource memory budget");
		}
		resource.second = p.memory_budget_bytes;
		resource.first = CheckedAdd(resource.first,
		    CheckedAdd(BasicDomainBytes(plan.Domains()[p.domain]), p.auxiliary_bytes));
		if (resource.first > resource.second)
		{
			throw std::invalid_argument(
			    "domain allocation exceeds aggregate resource memory budget");
		}
	}
}

bool ExchangeEpoch::operator==(const ExchangeEpoch& b) const
{
	return plan == b.plan && step == b.step && phase == b.phase;
}

bool MessageId::operator<(const MessageId& b) const
{
	return std::tie(epoch.plan, epoch.step, epoch.phase, channel, source, destination) <
	    std::tie(b.epoch.plan, b.epoch.step, b.epoch.phase, b.channel, b.source, b.destination);
}

bool MessageId::operator==(const MessageId& b) const
{
	return !(*this < b) && !(b < *this);
}

}
}

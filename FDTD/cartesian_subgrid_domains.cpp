/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesian_subgrid_domains.h"
#include "cartesian_subgrid.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace
{
using Index = std::uint32_t;
using Map = std::unordered_map<Index, Index>;
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}
void Rows(const std::vector<Index> &owned, const std::vector<Index> &starts,
	const std::vector<Index> &columns, const std::vector<float> &high, const std::vector<float> &low,
	const Map &input, std::vector<Index> &outStart, std::vector<Index> &outColumn,
	std::vector<float> &outHigh, std::vector<float> &outLow)
{
	outStart.assign(1, 0);
	for (Index row : owned)
	{
		for (Index k = starts[row]; k < starts[row + 1]; ++k)
		{
			// Mapping changes indices, never accumulation order or coefficients.
			outColumn.push_back(input.at(columns[k]));
			outHigh.push_back(high[k]);
			outLow.push_back(low[k]);
		}
		outStart.push_back(static_cast<Index>(outColumn.size()));
	}
}
}

CartesianSubgridDomainPlan::CartesianSubgridDomainPlan(const CartesianSubgrid &grid,
	unsigned int axis, const std::vector<double> &cuts)
{
	Require(axis < 3 && cuts.size() < 64, "CartesianSubgrid: invalid domain cuts");
	const auto &lines = grid.BaseLines()[axis];
	double previous = lines.front();
	for (double cut : cuts)
	{
		Require(std::isfinite(cut) && cut > previous && cut < lines.back(),
			"CartesianSubgrid: domain cuts must be finite, increasing and internal");
		previous = cut;
	}
	const auto data = grid.DeviceData();
	Require(data.rootConductance.empty() && data.portEdge.empty() && data.responseEdge.empty(),
		"CartesianSubgrid: distributed resistor ports are not implemented");
	m_StableStep = data.stableTimestep;
	m_FactorTails = data.factorTails;
	m_StateTails = data.stateTails;
	m_ElectricCount = data.electric.size();
	m_MagneticCount = data.magnetic.size();
	Require(m_MagneticCount == 6 * grid.CellCount(), "CartesianSubgrid: unexpected magnetic quadrature layout");
	const auto owner = [&](double coordinate)
	{
		return static_cast<Index>(std::upper_bound(cuts.begin(), cuts.end(), coordinate) - cuts.begin());
	};
	std::vector<Index> eOwner(m_ElectricCount), hOwner(m_MagneticCount);
	for (Index edge = 0; edge < m_ElectricCount; ++edge)
	{
		eOwner[edge] = owner(grid.EdgePosition(edge)[axis]);
	}
	for (Index leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const Index worker = owner(grid.LeafStart(leaf)[axis] + .5 * grid.LeafLengths(leaf)[axis]);
		std::fill(hOwner.begin() + 6 * leaf, hOwner.begin() + 6 * leaf + 6, worker);
	}
	for (Index block = 0; block < data.bandwidth.size(); ++block)
	{
		const Index first = data.blockStart[block], stop = data.blockStart[block + 1];
		double centre = 0;
		for (Index k = first; k < stop; ++k)
		{
			centre += grid.EdgePosition(data.blockEdge[k])[axis];
		}
		const Index worker = owner(centre / (stop - first));
		for (Index k = first; k < stop; ++k)
		{
			const Index edge = data.blockEdge[k];
			m_ReassignedElectric += eOwner[edge] != worker;
			eOwner[edge] = worker;
		}
	}
	m_Domains.resize(cuts.size() + 1);
	std::vector<Map> eMaps(m_Domains.size()), hMaps(m_Domains.size());
	for (Index worker = 0; worker < m_Domains.size(); ++worker)
	{
		auto &domain = m_Domains[worker];
		auto &local = domain.fields;
		local.stableTimestep = m_StableStep;
		for (Index i = 0; i < m_ElectricCount; ++i)
		{
			if (eOwner[i] == worker)
			{
				domain.globalElectric.push_back(i);
			}
		}
		for (Index i = 0; i < m_MagneticCount; ++i)
		{
			if (hOwner[i] == worker)
			{
				domain.globalMagnetic.push_back(i);
			}
		}
		domain.ownedElectric = domain.globalElectric.size();
		domain.ownedMagnetic = domain.globalMagnetic.size();
		Require(domain.ownedElectric && domain.ownedMagnetic, "CartesianSubgrid: empty domain is unsupported");
		const auto ownedE = domain.globalElectric, ownedH = domain.globalMagnetic;
		std::set<Index> ghostE, ghostH;
		for (Index row : ownedE)
		{
			for (Index k = data.transposeStart[row]; k < data.transposeStart[row + 1]; ++k)
			{
				const Index column = data.transposeColumn[k];
				if (hOwner[column] != worker)
				{
					ghostH.insert(column);
				}
			}
		}
		for (Index row : ownedH)
		{
			for (Index k = data.curlStart[row]; k < data.curlStart[row + 1]; ++k)
			{
				const Index column = data.curlColumn[k];
				if (eOwner[column] != worker)
				{
					ghostE.insert(column);
				}
			}
		}
		domain.globalElectric.insert(domain.globalElectric.end(), ghostE.begin(), ghostE.end());
		domain.globalMagnetic.insert(domain.globalMagnetic.end(), ghostH.begin(), ghostH.end());
		for (Index global : domain.globalElectric)
		{
			eMaps[worker].emplace(global, static_cast<Index>(local.electric.size()));
			local.electric.push_back(data.electric[global]);
		}
		for (Index global : domain.globalMagnetic)
		{
			hMaps[worker].emplace(global, static_cast<Index>(local.magnetic.size()));
			local.magnetic.push_back(data.magnetic[global]);
		}
		Rows(ownedE, data.transposeStart, data.transposeColumn, data.transposeWeight, data.transposeLow,
			hMaps[worker], local.transposeStart, local.transposeColumn, local.transposeWeight, local.transposeLow);
		Rows(ownedH, data.curlStart, data.curlColumn, data.curlWeight, data.curlLow,
			eMaps[worker], local.curlStart, local.curlColumn, local.curlWeight, local.curlLow);
		for (Index block = 0; block < data.bandwidth.size(); ++block)
		{
			const Index first = data.blockStart[block], stop = data.blockStart[block + 1];
			if (eOwner[data.blockEdge[first]] != worker)
			{
				continue;
			}
			for (Index k = first; k < stop; ++k)
			{
				Require(eOwner[data.blockEdge[k]] == worker, "CartesianSubgrid: split interface mass component");
				local.blockEdge.push_back(eMaps[worker].at(data.blockEdge[k]));
			}
			local.blockStart.push_back(static_cast<Index>(local.blockEdge.size()));
			local.bandwidth.push_back(data.bandwidth[block]);
			local.lower.insert(local.lower.end(), data.lower.begin() + data.lowerStart[block],
				data.lower.begin() + data.lowerStart[block + 1]);
			local.lowerStart.push_back(static_cast<Index>(local.lower.size()));
		}
		for (Index source = 0; source < data.sourceEdge.size(); ++source)
		{
			if (eOwner[data.sourceEdge[source]] == worker)
			{
				local.sourceEdge.push_back(eMaps[worker].at(data.sourceEdge[source]));
				local.sourceWeight.push_back(data.sourceWeight[source]);
				local.sourceOffStart.push_back(0);
			}
		}
	}
	for (bool electric : {true, false})
	{
		std::vector<Map> sends(m_Domains.size());
		auto &transfers = electric ? m_ElectricTransfers : m_MagneticTransfers;
		const auto &owners = electric ? eOwner : hOwner;
		const auto &maps = electric ? eMaps : hMaps;
		for (Index receiver = 0; receiver < m_Domains.size(); ++receiver)
		{
			const auto &domain = m_Domains[receiver];
			const auto &global = electric ? domain.globalElectric : domain.globalMagnetic;
			const std::size_t owned = electric ? domain.ownedElectric : domain.ownedMagnetic;
			for (Index k = static_cast<Index>(owned); k < global.size(); ++k)
			{
				const Index sender = owners[global[k]], source = maps[sender].at(global[k]);
				auto &indices = electric ? m_Domains[sender].sendElectric : m_Domains[sender].sendMagnetic;
				const auto inserted = sends[sender].emplace(source, static_cast<Index>(indices.size()));
				if (inserted.second)
				{
					indices.push_back(source);
				}
				transfers.push_back({sender, inserted.first->second, receiver, k - static_cast<Index>(owned)});
			}
		}
	}
}

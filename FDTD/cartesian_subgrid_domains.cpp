/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesian_subgrid_domains.h"
#include "cartesian_subgrid.h"
#include "subgrid_domain_map.h"
#include "extensions/engine_ext_cartesiansubgrid_conductivity.h"
#include "extensions/engine_ext_cartesiansubgrid_pml.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace
{
using Index = std::uint32_t;
using Map = std::unordered_map<Index, Index>;
// Setup-only row reader. The numerical grid already owns the exact curl.
// Query local integer support for C^T without global FP32 or incoming CSR.
// Ascending face traversal retains the established transpose sum order.
class CurlRows
{
  public:
	CurlRows(const CartesianSubgrid &grid, const CartesianSubgridDeviceData &data)
		: m_Grid(grid), m_Data(data)
	{
		m_Row.reserve(32);
	}

	template <class Visitor> void ForEach(bool transpose, Index row, Visitor visitor) const
	{
		if (!m_Data.streamCurls)
		{
			const auto &start = transpose ? m_Data.transposeStart : m_Data.curlStart;
			const auto &column = transpose ? m_Data.transposeColumn : m_Data.curlColumn;
			const auto &weight = transpose ? m_Data.transposeWeight : m_Data.curlWeight;
			const auto &low = transpose ? m_Data.transposeLow : m_Data.curlLow;
			for (Index k = start.at(row); k < start.at(std::size_t(row) + 1); ++k)
			{
				visitor(column[k], weight[k], low[k]);
			}
			return;
		}
		if (!transpose)
		{
			m_Grid.CurlRow(row, m_Row);
			for (const auto &term : m_Row)
			{
				Emit(static_cast<Index>(term.first), term.second, visitor);
			}
			return;
		}
		m_Grid.CurlTransposeRow(row, m_Row);
		for (const auto &term : m_Row)
		{
			Emit(static_cast<Index>(term.first), term.second, visitor);
		}
	}

  private:
	static float Scalar(double value)
	{
		const float converted = static_cast<float>(value);
		if (!std::isfinite(value) || !std::isfinite(converted) ||
			(value != 0 && std::abs(converted) < std::numeric_limits<float>::min()))
		{
			throw std::invalid_argument("CartesianSubgrid: nonrepresentable streamed curl coefficient");
		}
		return converted;
	}
	template <class Visitor> static void Emit(Index column, double value, const Visitor &visitor)
	{
		const float high = Scalar(value), low = Scalar(value - double(high));
		visitor(column, high, low);
	}
	const CartesianSubgrid &m_Grid;
	const CartesianSubgridDeviceData &m_Data;
	// Private to one serial preparation object; never shared by worker threads.
	mutable std::vector<std::pair<std::size_t, double>> m_Row;
};
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}
template <class Lookup> void Rows(const std::vector<Index> &owned, const std::vector<Index> &starts,
	const std::vector<Index> &columns, const std::vector<float> &high, const std::vector<float> &low,
	const Lookup &input, std::vector<Index> &outStart, std::vector<Index> &outColumn,
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
	unsigned int axis, const std::vector<double> &cuts, CurlPreparation curls)
	: CartesianSubgridDomainPlan(grid, axis, cuts, {{curls == CurlPreparation::ExpandedOracle ? grid.DeviceData() : grid.DeviceSetupMetadata(), {}}, {}, curls != CurlPreparation::ExpandedOracle})
{
}

CartesianSubgridDomainPlan::CartesianSubgridDomainPlan(
	const Engine_Ext_CartesianSubgrid_Conductivity &advance, unsigned int axis,
	const std::vector<double> &cuts, CurlPreparation curls)
	: CartesianSubgridDomainPlan(advance.Grid(), axis, cuts, {advance.PrepareDomainDataImpl(nullptr, curls == CurlPreparation::ExpandedOracle, curls == CurlPreparation::ExpandedOracle), {}, curls != CurlPreparation::ExpandedOracle})
{
}

CartesianSubgridDomainPlan::CartesianSubgridDomainPlan(const Engine_Ext_CartesianSubgrid_PML &pml,
	unsigned int axis, const std::vector<double> &cuts, CurlPreparation curls)
	: CartesianSubgridDomainPlan(pml.Grid(), axis, cuts, {{curls == CurlPreparation::ExpandedOracle ? pml.Grid().DeviceData() : pml.Grid().DeviceSetupMetadata(), {}},
		pml.PrepareDomainData(), curls != CurlPreparation::ExpandedOracle})
{
}

CartesianSubgridDomainPlan::CartesianSubgridDomainPlan(
	const Engine_Ext_CartesianSubgrid_Conductivity &advance, const Engine_Ext_CartesianSubgrid_PML &pml,
	unsigned int axis, const std::vector<double> &cuts, CurlPreparation curls)
	: CartesianSubgridDomainPlan(advance.Grid(), axis, cuts,
		{advance.PrepareDomainDataImpl(&pml, curls == CurlPreparation::ExpandedOracle, curls == CurlPreparation::ExpandedOracle), pml.PrepareDomainData(&advance), curls != CurlPreparation::ExpandedOracle})
{
}

CartesianSubgridDomainPlan::CartesianSubgridDomainPlan(const CartesianSubgrid &grid,
	unsigned int axis, const std::vector<double> &cuts, Prepared prepared)
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
	const auto &data = prepared.advance.fields;
	const CurlRows curls(grid, data);
	Require(!data.rootConductance.empty() || data.ledger == std::vector<float>({0, 0, 0, 0}),
		"CartesianSubgrid: a non-port accumulated ledger has no qualified ownership");
	m_StableStep = data.stableTimestep;
	m_PortStep = data.rootConductance.empty() ? 0 : data.portTimestep;
	m_RecurrenceStep = prepared.advance.loss ? prepared.advance.loss->timestep : 0;
	Require(!prepared.advance.loss || (m_RecurrenceStep > 0 && m_RecurrenceStep <= m_StableStep),
		"CartesianSubgrid: invalid domain recurrence timestep");
	m_PMLStep = prepared.pml ? prepared.pml->timestep : 0;
	Require(!prepared.pml || (m_PMLStep > 0 && m_PMLStep <= m_StableStep &&
		(m_RecurrenceStep == 0 || m_PMLStep == m_RecurrenceStep) &&
		(m_PortStep == 0 || m_PMLStep == m_PortStep)), "CartesianSubgrid: invalid domain PML timestep");
	m_FactorTails = data.factorTails;
	m_StateTails = data.stateTails;
	Require(!prepared.borrowedInitialState ||
		(data.streamCurls && data.electric.empty() && data.magnetic.empty()),
		"CartesianSubgrid: borrowed initial state requires checked metadata only");
	m_ElectricCount = prepared.borrowedInitialState ? grid.EdgeCount() : data.electric.size();
	m_MagneticCount = prepared.borrowedInitialState ? grid.FaceCount() : data.magnetic.size();
	Require(m_MagneticCount == 6 * grid.CellCount(), "CartesianSubgrid: unexpected magnetic quadrature layout");
	const auto owner = [&](double coordinate)
	{
		return static_cast<Index>(std::upper_bound(cuts.begin(), cuts.end(), coordinate) - cuts.begin());
	};
	CartesianSubgridDomainOwners eOwner(data.streamCurls), hOwner(data.streamCurls);
	eOwner.reserve(m_ElectricCount);
	hOwner.reserve(m_MagneticCount);
	for (Index edge = 0; edge < m_ElectricCount; ++edge)
	{
		eOwner.push_back(owner(grid.EdgePosition(edge)[axis]));
	}
	for (Index leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const Index worker = owner(grid.LeafStart(leaf)[axis] + .5 * grid.LeafLengths(leaf)[axis]);
		for (unsigned int face = 0; face < 6; ++face)
		{
			hOwner.push_back(worker);
		}
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
			eOwner.Set(edge, worker);
		}
	}
	// Preserve the complete small Woodbury solve, its ordering and responses.
	// Co-locate all ports and every touched mass component; only those electric
	// states move. Independent port-cluster splitting is a later optimization.
	Index portOwner = 0;
	if (!data.rootConductance.empty())
	{
		std::vector<bool> touched(m_ElectricCount, false);
		for (Index edge : data.portEdge)
		{
			touched[edge] = true;
		}
		for (Index edge : data.responseEdge)
		{
			touched[edge] = true;
		}
		for (std::size_t block = 0; block < data.bandwidth.size(); ++block)
		{
			const auto first = data.blockEdge.begin() + data.blockStart[block];
			const auto stop = data.blockEdge.begin() + data.blockStart[block + 1];
			if (std::any_of(first, stop, [&](Index edge) { return touched[edge]; }))
			{
				for (auto edge = first; edge != stop; ++edge)
				{
					touched[*edge] = true;
				}
			}
		}
		double centre = 0;
		std::size_t count = 0;
		for (Index edge = 0; edge < m_ElectricCount; ++edge)
		{
			if (touched[edge])
			{
				centre += grid.EdgePosition(edge)[axis];
				++count;
			}
		}
		Require(count != 0, "CartesianSubgrid: ports have no owned electric support");
		portOwner = owner(centre / count);
		for (Index edge = 0; edge < m_ElectricCount; ++edge)
		{
			if (touched[edge])
			{
				eOwner.Set(edge, portOwner);
			}
		}
	}
	for (Index edge = 0; edge < m_ElectricCount; ++edge)
	{
		m_ReassignedElectric += eOwner[edge] != owner(grid.EdgePosition(edge)[axis]);
	}
	m_Domains.resize(cuts.size() + 1);
	// One direct local address for every globally owned state. Only ghosts
	// need a per-worker hash entry. The expanded oracle retains full maps.
	std::vector<Index> eLocal, hLocal;
	if (data.streamCurls)
	{
		eLocal.assign(m_ElectricCount, CartesianSubgridDomainMap::absent);
		hLocal.assign(m_MagneticCount, CartesianSubgridDomainMap::absent);
	}
	std::vector<CartesianSubgridDomainMap> eMaps, hMaps;
	eMaps.reserve(m_Domains.size());
	hMaps.reserve(m_Domains.size());
	for (Index worker = 0; worker < m_Domains.size(); ++worker)
	{
		eMaps.emplace_back(eOwner, eLocal, worker, data.streamCurls);
		hMaps.emplace_back(hOwner, hLocal, worker, data.streamCurls);
	}
	for (Index worker = 0; worker < m_Domains.size(); ++worker)
	{
		auto &domain = m_Domains[worker];
		auto &local = domain.fields;
		local.stableTimestep = m_StableStep;
		local.trackSourceWork = data.trackSourceWork || !data.rootConductance.empty();
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
			curls.ForEach(true, row, [&](Index column, float, float)
			{
				if (hOwner[column] != worker)
				{
					ghostH.insert(column);
				}
			});
		}
		for (Index row : ownedH)
		{
			curls.ForEach(false, row, [&](Index column, float, float)
			{
				if (eOwner[column] != worker)
				{
					ghostE.insert(column);
				}
			});
		}
		if (prepared.pml)
		{
			const auto dependencies = [&](const CartesianSubgridPMLData::Rows &rows,
				const CartesianSubgridDomainOwners &targetOwners, const CartesianSubgridDomainOwners &sourceOwners,
				std::set<Index> &ghosts)
			{
				for (Index row = 0; row < rows.target.size(); ++row)
				{
					if (targetOwners.at(rows.target[row]) != worker)
					{
						continue;
					}
					for (Index k = rows.start.at(row); k < rows.start.at(row + 1); ++k)
					{
						const Index source = rows.source.at(k);
						if (sourceOwners.at(source) != worker)
						{
							ghosts.insert(source);
						}
					}
				}
			};
			dependencies(prepared.pml->electric, eOwner, hOwner, ghostH);
			dependencies(prepared.pml->magnetic, hOwner, eOwner, ghostE);
		}
		domain.globalElectric.insert(domain.globalElectric.end(), ghostE.begin(), ghostE.end());
		domain.globalMagnetic.insert(domain.globalMagnetic.end(), ghostH.begin(), ghostH.end());
		for (Index global : domain.globalElectric)
		{
			eMaps[worker].emplace(global, static_cast<Index>(local.electric.size()));
			local.electric.push_back(prepared.borrowedInitialState ?
				CartesianSubgridFloat::StateValue(grid.m_Electric.at(global)) : data.electric[global]);
		}
		for (Index global : domain.globalMagnetic)
		{
			hMaps[worker].emplace(global, static_cast<Index>(local.magnetic.size()));
			local.magnetic.push_back(prepared.borrowedInitialState ?
				CartesianSubgridFloat::StateValue(grid.m_Magnetic.at(global)) : data.magnetic[global]);
		}
		const auto rows = [&](bool transpose, const std::vector<Index> &owned, const CartesianSubgridDomainMap &input,
			std::vector<Index> &starts, std::vector<Index> &columns,
			std::vector<float> &weight, std::vector<float> &low)
		{
			using Curl = CartesianSubgridDeviceCurl;
			std::optional<Curl::Builder> stream;
			std::vector<Curl::Term> terms;
			if (data.streamCurls)
			{
				stream.emplace(owned.size(), input.size(), true);
				terms.reserve(8);
			}
			else
			{
				starts.assign(1, 0);
			}
			for (Index row : owned)
			{
				terms.clear();
				curls.ForEach(transpose, row, [&](Index column, float high, float residual)
				{
					const Index localColumn = input.at(column);
					if (stream)
					{
						terms.push_back({localColumn, high, residual});
					}
					else
					{
						columns.push_back(localColumn);
						weight.push_back(high);
						low.push_back(residual);
					}
				});
				if (stream)
				{
					stream->Append(terms);
				}
				else
				{
					starts.push_back(static_cast<Index>(columns.size()));
				}
			}
			if (stream)
			{
				if (!domain.operators)
				{
					domain.operators.emplace();
				}
				auto &result = transpose ? domain.operators->transpose : domain.operators->curl;
				result = stream->Finish();
				local.streamCurls = true;
			}
		};
		rows(true, ownedE, hMaps[worker], local.transposeStart, local.transposeColumn,
			local.transposeWeight, local.transposeLow);
		rows(false, ownedH, eMaps[worker], local.curlStart, local.curlColumn, local.curlWeight, local.curlLow);
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
				for (Index k = data.sourceOffStart[source]; k < data.sourceOffStart[source + 1]; ++k)
				{
					Require(eOwner[data.sourceOffEdge[k]] == worker,
						"CartesianSubgrid: source work crosses a split mass component");
					local.sourceOffEdge.push_back(eMaps[worker].at(data.sourceOffEdge[k]));
					local.sourceOffWeight.push_back(data.sourceOffWeight[k]);
				}
				local.sourceOffStart.push_back(static_cast<Index>(local.sourceOffEdge.size()));
			}
		}
		if (!data.rootConductance.empty() && worker == portOwner)
		{
			local.portTimestep = data.portTimestep;
			local.rootConductance = data.rootConductance;
			local.portCurrent = data.portCurrent;
			local.portLower = data.portLower;
			local.portStart = data.portStart;
			local.responseStart = data.responseStart;
			local.portWeight = data.portWeight;
			local.responseWeight = data.responseWeight;
			for (Index edge : data.portEdge)
			{
				Require(eOwner[edge] == worker, "CartesianSubgrid: split port voltage support");
				local.portEdge.push_back(eMaps[worker].at(edge));
			}
			for (Index edge : data.responseEdge)
			{
				Require(eOwner[edge] == worker, "CartesianSubgrid: split port response support");
				local.responseEdge.push_back(eMaps[worker].at(edge));
			}
			// Attribute the existing accumulated work once, never per worker.
			local.ledger = data.ledger;
		}
		if (prepared.pml)
		{
			domain.pml.emplace();
			domain.pml->timestep = prepared.pml->timestep;
			const auto select = [&](const CartesianSubgridPMLData::Rows &input,
				const CartesianSubgridDomainOwners &owners, const auto &targets, const auto &sources,
				CartesianSubgridPMLData::Rows &output)
			{
				for (Index group = 0; group + 1 < input.group.size(); ++group)
				{
					const Index first = input.group[group], stop = input.group[group + 1];
					if (owners.at(input.target.at(first)) != worker)
					{
						continue;
					}
					for (Index row = first; row < stop; ++row)
					{
						Require(owners.at(input.target[row]) == worker, "CartesianSubgrid: split PML target group");
						output.target.push_back(targets.at(input.target[row]));
						output.decay.push_back(input.decay[row]);
						output.change.push_back(input.change[row]);
						for (Index k = input.start[row]; k < input.start[row + 1]; ++k)
						{
							output.source.push_back(sources.at(input.source[k]));
							output.weight.push_back(input.weight[k]);
							output.low.push_back(input.low[k]);
						}
						output.start.push_back(static_cast<Index>(output.source.size()));
					}
					output.group.push_back(static_cast<Index>(output.target.size()));
				}
			};
			select(prepared.pml->electric, eOwner, eMaps[worker], hMaps[worker], domain.pml->electric);
			select(prepared.pml->magnetic, hOwner, hMaps[worker], eMaps[worker], domain.pml->magnetic);
		}
		if (prepared.advance.loss)
		{
			const auto &source = *prepared.advance.loss;
			domain.loss.emplace();
			auto &loss = *domain.loss;
			loss.timestep = source.timestep;
			for (Index row : ownedH)
			{
				loss.magneticRate.push_back(source.magneticRate.at(row));
			}
			const auto ownedRows = [&](const CartesianSubgridConductivityData::Rows &input,
				const std::vector<Index> &rows, const auto &columns, CartesianSubgridConductivityData::Rows &output)
			{
				Rows(rows, input.start, input.edge, input.weight, input.low,
					columns, output.start, output.edge, output.weight, output.low);
			};
			// Every conductance and surface trace must stay within one complete
			// electric mass component. Ghost fields are for curls, never histories.
			const auto ownedMap = eMaps[worker].Owned();
			ownedRows(source.conductance, ownedE, ownedMap, loss.conductance);
			const auto select = [&](const CartesianSubgridConductivityData::Rows &rows)
			{
				std::vector<Index> selected;
				for (Index row = 0; row + 1 < rows.start.size(); ++row)
				{
					const Index first = rows.start[row], stop = rows.start[row + 1];
					Require(first < stop, "CartesianSubgrid: empty constitutive trace");
					const Index assigned = eOwner.at(rows.edge.at(first));
					for (Index k = first; k < stop; ++k)
					{
						Require(eOwner.at(rows.edge.at(k)) == assigned,
							"CartesianSubgrid: constitutive trace crosses ownership");
					}
					if (assigned == worker)
					{
						selected.push_back(row);
					}
				}
				return selected;
			};
			ownedRows(source.trace, select(source.trace), ownedMap, loss.trace);
			const auto poles = select(source.surfaceTrace);
			if (!poles.empty())
			{
				Map poleMap;
				for (Index k = 0; k < poles.size(); ++k)
				{
					poleMap.emplace(poles[k], k);
					loss.surfaceRate.push_back(source.surfaceRate.at(poles[k]));
				}
				ownedRows(source.surfaceTrace, poles, ownedMap, loss.surfaceTrace);
				ownedRows(source.surfaceHistory, ownedE, poleMap, loss.surfaceHistory);
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

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "operator_ext_cartesiansubgrid_mur.h"
#include "CSPrimitives.h"
#include "CSPropExcitation.h"
#include "ContinuousStructure.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "engine_ext_cartesiansubgrid_mur.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
void Require(bool valid, const char *message)
{
	if (!valid)
	{
		throw std::invalid_argument(message);
	}
}
} // namespace

Operator_Ext_CartesianSubgrid_Mur::Operator_Ext_CartesianSubgrid_Mur(Operator *op, unsigned int face)
	: Operator_Extension(op), m_Compact(dynamic_cast<Operator_CartesianSubgrid *>(op)), m_Face(face)
{
	Require(m_Compact && face < 6, "CartesianSubgrid: invalid Mur operator or face");
}

void Operator_Ext_CartesianSubgrid_Mur::SetPhaseVelocity(double velocity)
{
	Require(velocity > 0 && std::isfinite(velocity),
			"CartesianSubgrid: finite positive Mur velocity required");
	m_Velocity = velocity;
}

bool Operator_Ext_CartesianSubgrid_Mur::BuildExtension()
{
	auto &grid = m_Compact->Grid();
	const auto &axes = grid.BaseLines();
	const unsigned int normal = m_Face / 2;
	const bool high = m_Face % 2;
	const auto &axis = axes[normal];
	const std::size_t end = axis.size() - 1;
	const double boundary = high ? axis.back() : axis.front();
	const double inside = high ? axis[end - 1] : axis[1];
	const double collar = high ? axis[end - 2] : axis[2];
	const double tolerance = 64 * std::numeric_limits<double>::epsilon() *
							 std::max(std::abs(boundary), std::abs(boundary - inside));
	for (const auto &region : grid.Regions())
	{
		Require(high ? region.box.stop[normal] <= end - 2 : region.box.start[normal] >= 2,
				"CartesianSubgrid: Mur requires two unrefined base cells at its face");
	}
	CartesianSubgrid::Medium medium{};
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
		if (high ? lo[normal] + length[normal] <= collar + tolerance : lo[normal] >= collar - tolerance)
		{
			continue;
		}
		const auto value = grid.LeafMedium(leaf);
		if (medium.epsilon == 0)
		{
			medium = value;
		}
		Require(std::abs(value.epsilon / medium.epsilon - 1) < 1e-12 &&
					std::abs(value.permeability / medium.permeability - 1) < 1e-12,
				"CartesianSubgrid: Mur collar must be homogeneous and lossless");
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int a = 0; a < 2; ++a)
			{
				for (unsigned int b = 0; b < 2; ++b)
				{
					Require(!grid.LeafEdgeIsPEC(leaf, n, a, b),
							"CartesianSubgrid: internal PEC touches the Mur collar");
				}
			}
		}
	}
	for (auto *property : m_Op->GetGeometryCSX()->GetPropertyByType(CSProperties::ANY))
	{
		if (property->GetType() != CSProperties::EXCITATION &&
			property->GetType() != CSProperties::LUMPED_ELEMENT)
		{
			continue;
		}
		if (auto *source = dynamic_cast<CSPropExcitation *>(property))
		{
			if (!source->GetEnabled())
			{
				continue;
			}
		}
		for (auto *primitive : property->GetAllPrimitives())
		{
			double box[6];
			primitive->GetBoundBox(box);
			double start[3] = {box[0], box[2], box[4]}, stop[3] = {box[1], box[3], box[5]};
			const auto a = m_Compact->DrawingPoint(start), b = m_Compact->DrawingPoint(stop);
			Require(high ? b[normal] < collar - tolerance : a[normal] > collar + tolerance,
					"CartesianSubgrid: sources and resistors must stay outside the Mur collar");
		}
	}
	Require(medium.epsilon > 0, "CartesianSubgrid: empty Mur collar");
	const double velocity = m_Velocity > 0 ? m_Velocity : 1 / std::sqrt(medium.epsilon * medium.permeability);
	const double travel = velocity * m_Op->GetTimestep(), distance = std::abs(boundary - inside);
	m_Coefficient = (travel - distance) / (travel + distance);
	Require(std::isfinite(m_Coefficient) && std::abs(m_Coefficient) < 1,
			"CartesianSubgrid: invalid Mur coefficient");
	m_Nodes.clear();
	for (std::size_t edge = 0; edge < grid.EdgeCount(); ++edge)
	{
		const unsigned int n = grid.EdgeDirection(edge);
		const auto centre = grid.EdgePosition(edge);
		if (n == normal || std::abs(centre[normal] - boundary) > tolerance)
		{
			continue;
		}
		const auto &tangent = axes[n];
		const auto after = std::upper_bound(tangent.begin(), tangent.end(), centre[n]);
		Require(after != tangent.begin() && after != tangent.end(),
				"CartesianSubgrid: invalid Mur edge location");
		auto start = centre, stop = centre;
		start[n] = *(after - 1);
		stop[n] = *after;
		const auto trace = grid.VoltagePath(start, stop);
		start[normal] = stop[normal] = inside;
		const auto neighbour = grid.VoltagePath(start, stop);
		Require(trace.weights.size() == 1 && trace.weights[0].first == edge && neighbour.weights.size() == 1,
				"CartesianSubgrid: Mur requires independent boundary and neighbour edges");
		const double weight = trace.weights[0].second;
		m_Nodes.push_back({trace, neighbour, edge, stop[n] - start[n], 1 / (weight * weight)});
	}
	Require(!m_Nodes.empty(), "CartesianSubgrid: no independent Mur boundary traces");
	return true;
}

Engine_Extension *Operator_Ext_CartesianSubgrid_Mur::CreateEngineExtention()
{
	m_Eng_Ext = new Engine_Ext_CartesianSubgrid_Mur(this);
	return m_Eng_Ext;
}

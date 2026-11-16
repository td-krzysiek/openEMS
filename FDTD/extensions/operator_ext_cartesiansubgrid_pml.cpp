/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "operator_ext_cartesiansubgrid_pml.h"
#include "CSPrimitives.h"
#include "CSPropExcitation.h"
#include "ContinuousStructure.h"
#include "FDTD/operator_cartesiansubgrid.h"
#include "engine_ext_cartesiansubgrid_pml.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace
{
void Require(bool valid, const char *message)
{
	if (!valid)
		throw std::invalid_argument(message);
}
} // namespace

Operator_Ext_CartesianSubgrid_PML::Operator_Ext_CartesianSubgrid_PML(Operator *op,
																	 const std::array<unsigned int, 6> &cells)
	: Operator_Extension(op), m_Compact(dynamic_cast<Operator_CartesianSubgrid *>(op)), m_Cells(cells)
{
	Require(m_Compact && std::any_of(cells.begin(), cells.end(), [](unsigned int n) { return n != 0; }),
			"CartesianSubgrid: PML requires a compact operator and at least one face");
}

bool Operator_Ext_CartesianSubgrid_PML::BuildExtension()
{
	m_Prepared = false;
	m_Footprint = {};
	auto &g = m_Compact->Grid();
	const auto &axes = g.BaseLines();
	std::array<double, 6> inner{}, guard{}, width{};
	const double speed = 1 / std::sqrt(g.m_Epsilon * g.m_Mu);
	const double dt = m_Op->GetTimestep();
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (!m_Cells[2 * n] && !m_Cells[2 * n + 1])
			continue;
		const unsigned int end = unsigned(axes[n].size() - 1);
		Require(std::size_t(m_Cells[2 * n]) + m_Cells[2 * n + 1] + 2 <= end,
				"CartesianSubgrid: opposing PML layers need two intervening base cells");
		for (unsigned int side = 0; side < 2; ++side)
		{
			const unsigned int face = 2 * n + side, count = m_Cells[face];
			if (!count)
				continue;
			Require(g.m_Boundaries[face] == CartesianSubgrid::BoundaryType::PEC,
					"CartesianSubgrid: PML requires a PEC termination at its outer face");
			Require(count >= 4, "CartesianSubgrid: PML requires at least four base cells per face");
			inner[face] = axes[n][side ? end - count : count];
			guard[face] = axes[n][side ? end - count - 1 : count + 1];
			width[face] = side ? axes[n].back() - inner[face] : inner[face] - axes[n].front();
			for (const auto &region : g.Regions())
				Require(side ? region.box.stop[n] <= end - count - 1 : region.box.start[n] >= count + 1,
						"CartesianSubgrid: refinement must stay outside the PML plus one-cell guard");
		}
	}
	auto collar = [&](const CartesianSubgrid::Point &lo, const CartesianSubgrid::Point &hi)
	{
		for (unsigned int face = 0; face < 6; ++face)
		{
			if (!m_Cells[face])
				continue;
			const unsigned int n = face / 2;
			const double tol =
				64 * std::numeric_limits<double>::epsilon() *
				std::max({std::abs(guard[face]), width[face], axes[n].back() - axes[n].front()});
			if (face % 2 ? hi[n] > guard[face] + tol : lo[n] < guard[face] - tol)
				return true;
		}
		return false;
	};
	for (std::size_t leaf = 0; leaf < g.CellCount(); ++leaf)
	{
		const auto lo = g.LeafStart(leaf), length = g.LeafLengths(leaf);
		auto hi = lo;
		for (unsigned int n = 0; n < 3; ++n)
			hi[n] += length[n];
		if (!collar(lo, hi))
			continue;
		const auto material = g.LeafMedium(leaf);
		Require(std::abs(material.epsilon / g.m_Epsilon - 1) < 1e-12 &&
					std::abs(material.permeability / g.m_Mu - 1) < 1e-12,
				"CartesianSubgrid: PML collar must match the homogeneous lossless background");
		for (std::size_t face = 6 * leaf; face < 6 * (leaf + 1); ++face)
		{
			m_Footprint.magnetic.push_back(face);
		}
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int a = 0; a < 2; ++a)
			{
				for (unsigned int b = 0; b < 2; ++b)
				{
					const auto trace = g.Restrict(g.LocalEdge(g.m_Leaves[leaf], n, a, b));
					for (const auto &term : trace.weights)
					{
						m_Footprint.electric.push_back(term.first);
					}
				}
			}
		}
		for (unsigned int n = 0; n < 3; ++n)
			for (unsigned int a = 0; a < 2; ++a)
				for (unsigned int b = 0; b < 2; ++b)
					Require(!g.LeafEdgeIsPEC(leaf, n, a, b),
							"CartesianSubgrid: internal PEC touches the PML collar");
	}
	for (auto *property : m_Op->GetGeometryCSX()->GetPropertyByType(CSProperties::ANY))
	{
		if (property->GetType() != CSProperties::EXCITATION &&
			property->GetType() != CSProperties::LUMPED_ELEMENT)
			continue;
		if (auto *source = dynamic_cast<CSPropExcitation *>(property))
			if (!source->GetEnabled())
				continue;
		for (auto *primitive : property->GetAllPrimitives())
		{
			double box[6];
			primitive->GetBoundBox(box);
			double a[3] = {box[0], box[2], box[4]}, b[3] = {box[1], box[3], box[5]};
			const auto lo = m_Compact->DrawingPoint(a), hi = m_Compact->DrawingPoint(b);
			for (unsigned int face = 0; face < 6; ++face)
				if (m_Cells[face])
					Require(face % 2 ? hi[face / 2] < guard[face] : lo[face / 2] > guard[face],
							"CartesianSubgrid: sources and resistors must stay outside the PML collar");
		}
	}
	auto rate = [&](unsigned int axis, double coordinate)
	{
		double value = 0;
		for (unsigned int side = 0; side < 2; ++side)
		{
			const unsigned int face = 2 * axis + side;
			if (!m_Cells[face])
				continue;
			const double depth = (side ? coordinate - inner[face] : inner[face] - coordinate) / width[face];
			// Cubic profile. R=1e-8 is a continuous round-trip amplitude
			// design parameter, not a claim about discrete measured reflection.
			if (depth > 0)
				value += -2 * speed * std::log(1e-8) / width[face] * std::pow(depth, 3);
		}
		return value;
	};
	using Key = std::pair<std::size_t, unsigned int>;
	std::map<Key, Derivative> electric, magnetic;
	auto add = [&](std::map<Key, Derivative> &nodes, std::size_t target, unsigned int axis, double coordinate,
				   std::size_t source, double coefficient)
	{
		const double damping = rate(axis, coordinate);
		if (damping <= 0)
			return;
		auto result = nodes.emplace(
			Key{target, axis}, Derivative{target, std::exp(-damping * dt), std::expm1(-damping * dt), {}});
		result.first->second.terms.push_back({source, coefficient});
	};
	for (std::size_t row = 0; row < g.FaceCount(); ++row)
	{
		const std::size_t leaf = row / 6;
		const unsigned int normal = unsigned(row % 6) / 2, side = unsigned(row % 2);
		auto point = g.LeafStart(leaf);
		const auto length = g.LeafLengths(leaf);
		for (unsigned int n = 0; n < 3; ++n)
			point[n] += length[n] * (n == normal ? double(side) : .5);
		g.m_CurlOperator.ForEach(row, [&](std::size_t edge, double coefficient)
		{
			const unsigned int tangent = g.EdgeDirection(edge);
			Require(normal != tangent, "CartesianSubgrid: invalid PML curl orientation");
			const unsigned int axis = 3 - normal - tangent;
			add(magnetic, row, axis, point[axis], edge, coefficient);
			add(electric, edge, axis, g.EdgePosition(edge)[axis], row, coefficient);
		});
	}
	// Collar corrections cannot bypass the non-diagonal interface mass.
	for (const auto &block : g.m_Blocks)
		for (std::size_t edge : block.edges)
			for (unsigned int n = 0; n < 3; ++n)
				Require(!electric.count(Key{edge, n}),
						"CartesianSubgrid: PML correction touches coupled interface mass");
	m_Electric.clear();
	m_Magnetic.clear();
	for (auto &node : electric)
		m_Electric.push_back(std::move(node.second));
	for (auto &node : magnetic)
		m_Magnetic.push_back(std::move(node.second));
	Require(!m_Electric.empty() && !m_Magnetic.empty(), "CartesianSubgrid: empty PML derivative state");
	for (const auto &node : m_Electric)
	{
		m_Footprint.electric.push_back(node.target);
	}
	for (const auto &node : m_Magnetic)
	{
		m_Footprint.magnetic.push_back(node.target);
	}
	for (auto *targets : {&m_Footprint.electric, &m_Footprint.magnetic})
	{
		std::sort(targets->begin(), targets->end());
		targets->erase(std::unique(targets->begin(), targets->end()), targets->end());
	}
	m_Prepared = true;
	return true;
}

const Operator_Ext_CartesianSubgrid_PML::Footprint &
Operator_Ext_CartesianSubgrid_PML::PreparedFootprint() const
{
	Require(m_Prepared, "CartesianSubgrid: PML footprint requested before successful preparation");
	return m_Footprint;
}

Engine_Extension *Operator_Ext_CartesianSubgrid_PML::CreateEngineExtention()
{
	m_Eng_Ext = new Engine_Ext_CartesianSubgrid_PML(this);
	return m_Eng_Ext;
}

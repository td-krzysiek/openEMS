/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesiansubgrid_conductor.h"
#include "surface_impedance.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <limits>
#include <tuple>
#include <stdexcept>
#include <utility>

namespace CartesianSubgridConductor
{
Prepared Prepare(const CartesianSubgrid &grid, const Sampler &sample, double maximumFrequency)
{
	if (!sample || !(maximumFrequency > 0) || !std::isfinite(maximumFrequency))
	{
		throw std::invalid_argument("CartesianSubgrid: smooth conductor needs a sampler and finite positive fit frequency");
	}
	Prepared result;
	struct Run
	{
		double lower, upper, conductivity;
	};
	using ColumnKey = std::tuple<unsigned int, double, double>;
	std::map<ColumnKey, std::vector<Run>> columns;
	std::vector<double> conductivity(grid.CellCount(), std::numeric_limits<double>::quiet_NaN());
	const auto readConductivity = [&](std::size_t leaf)
	{
		double &value = conductivity[leaf];
		if (std::isnan(value))
		{
			value = sample(leaf);
		}
		if (!(value > 0) || !std::isfinite(value))
		{
			throw std::invalid_argument("CartesianSubgrid: finite positive opaque-leaf conductivity required");
		}
		return value;
	};
	// A ray uses actual occupied cells, not the first containing primitive.
	// Cache each projected column so overlapping/touching box decompositions
	// give the same normal thickness without repeating CSXCAD queries per face.
	const auto normalRun = [&](std::size_t leaf, unsigned int axis) -> Run
	{
		const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
		auto point = lo;
		for (unsigned int n = 0; n < 3; ++n)
		{
			point[n] += .5 * length[n];
		}
		const ColumnKey key{axis, point[(axis + 1) % 3], point[(axis + 2) % 3]};
		auto found = columns.find(key);
		if (found == columns.end())
		{
			const auto &lines = grid.BaseLines()[axis];
			std::vector<Run> runs;
			double position = lines.front();
			bool previousOpaque = false;
			while (position < lines.back())
			{
				auto query = point;
				// FindLeaf already defines the positive-side tie at an exact face.
				// Do not introduce a global epsilon that can skip thin remote cells.
				query[axis] = position;
				const auto current = grid.LeafAt(query);
				const auto start = grid.LeafStart(current), size = grid.LeafLengths(current);
				const double end = start[axis] + size[axis];
				if (!(end > position))
				{
					throw std::invalid_argument("CartesianSubgrid: opaque normal-column walk made no progress");
				}
				if (grid.LeafIsOpaque(current))
				{
					const double sigma = readConductivity(current);
					if (previousOpaque)
					{
						if (sigma != runs.back().conductivity)
						{
							throw std::invalid_argument("CartesianSubgrid: touching heterogeneous surface conductors need a separate interface model");
						}
						runs.back().upper = end;
					}
					else
					{
						runs.push_back({start[axis], end, sigma});
					}
					previousOpaque = true;
				}
				else
				{
					previousOpaque = false;
				}
				position = end;
			}
			found = columns.emplace(key, std::move(runs)).first;
		}
		for (const auto &run : found->second)
		{
			if (point[axis] > run.lower && point[axis] < run.upper)
			{
				return run;
			}
		}
		throw std::invalid_argument("CartesianSubgrid: surface trace has no occupied normal-column run");
	};
	const auto corners = SurfaceImpedance::CornerFrequencies(maximumFrequency, 8);
	std::map<std::pair<double, double>, SurfaceImpedance::Fit> fits;
	for (const auto &face : grid.OpaqueBoundaryFaces())
	{
		const auto run = normalRun(face.leaf, face.normal);
		const double sigma = run.conductivity, thickness = run.upper - run.lower;
		if (!(sigma > 0) || !std::isfinite(sigma) || !(thickness > 0) || !std::isfinite(thickness))
		{
			throw std::invalid_argument("CartesianSubgrid: finite positive surface conductivity/body thickness required");
		}
		const auto key = std::make_pair(sigma, thickness);
		auto found = fits.find(key);
		if (found == fits.end())
		{
			auto fit = SurfaceImpedance::FitFace(sigma, thickness, 0, maximumFrequency, corners);
			if (fit.d_ref != 0 || fit.a.size() != corners.size() || !std::isfinite(fit.G) || fit.G < 0 ||
				!std::isfinite(fit.err_re) || !std::isfinite(fit.err_im) || fit.err_re < 0 || fit.err_im < 0)
			{
				throw std::invalid_argument("CartesianSubgrid: invalid or shifted smooth surface fit");
			}
			result.maximumRealFitError = std::max(result.maximumRealFitError, fit.err_re);
			result.maximumImaginaryFitError = std::max(result.maximumImaginaryFitError, fit.err_im);
			found = fits.emplace(key, std::move(fit)).first;
		}
		const auto &fit = found->second;
		const auto origin = grid.LeafStart(face.leaf), length = grid.LeafLengths(face.leaf);
		const double area = length[(face.normal + 1) % 3] * length[(face.normal + 2) % 3];
		for (unsigned int tangent = 1; tangent <= 2; ++tangent)
		{
			const unsigned int axis = (face.normal + tangent) % 3;
			const double weight = area / (2 * length[axis] * length[axis]);
			if (!(weight > 0) || !std::isfinite(weight))
			{
				throw std::invalid_argument("CartesianSubgrid: unrepresentable surface quadrature weight");
			}
			for (unsigned int transverse = 0; transverse < 2; ++transverse)
			{
				using Surface = Operator_Ext_CartesianSubgrid_Conductivity::SurfaceAdmittance;
				Surface surface;
				Surface::LocalEdge edge{face.leaf, axis, 0, 0};
				edge.first = (axis + 1) % 3 == face.normal ? face.side : transverse;
				edge.second = (axis + 2) % 3 == face.normal ? face.side : transverse;
				surface.localEdge = edge;
				surface.start = origin;
				surface.start[(axis + 1) % 3] += edge.first * length[(axis + 1) % 3];
				surface.start[(axis + 2) % 3] += edge.second * length[(axis + 2) % 3];
				surface.stop = surface.start;
				surface.stop[axis] += length[axis];
				if (grid.LeafEdgeVoltage(edge.leaf, edge.axis, edge.first, edge.second).weights.empty())
				{
					continue; // An independently imposed PEC contact/outer boundary.
				}
				surface.conductance = weight * fit.G;
				for (std::size_t pole = 0; pole < corners.size(); ++pole)
				{
					if (!std::isfinite(fit.a[pole]) || fit.a[pole] < 0)
					{
						throw std::invalid_argument("CartesianSubgrid: non-passive surface fit coefficient");
					}
					if (fit.a[pole] > 0)
					{
						const double resistance = 1 / (weight * fit.a[pole]);
						surface.poles.push_back({resistance, resistance / corners[pole]});
					}
				}
				result.surfaces.push_back(std::move(surface));
			}
		}
		++result.faces;
	}
	result.fits = fits.size();
	result.normalColumns = columns.size();
	return result;
}
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesian_subgrid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace
{
void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::invalid_argument(message);
	}
}
CartesianSubgrid::Axes UniformLines(CartesianSubgrid::Index cells, CartesianSubgrid::Point spacing)
{
	CartesianSubgrid::Axes lines;
	for (unsigned int n = 0; n < 3; ++n)
	{
		Require(cells[n] >= 4 && cells[n] <= 1024 && spacing[n] > 0 && std::isfinite(spacing[n]),
				"CartesianSubgrid: invalid base dimensions or spacing");
		for (unsigned int i = 0; i <= cells[n]; ++i)
		{
			lines[n].push_back(i * spacing[n]);
		}
	}
	return lines;
}

void CheckSubdivision(double lo, double hi, unsigned int ratio)
{
	const double width = hi - lo;
	// Sampling uses edge/face centres too. Reject unrepresentable local
	// positions rather than silently collapsing distinct cells or centres.
	double previous = lo;
	for (unsigned int i = 1; i <= 2 * ratio; ++i)
	{
		const double next = i == 2 * ratio ? hi : lo + width * (double(i) / (2 * ratio));
		Require(std::isfinite(next) && next > previous,
				"CartesianSubgrid: subdivision or sample centre is not representable at this coordinate scale");
		previous = next;
	}
}
} // namespace

bool CartesianSubgrid::Edge::operator<(const Edge &other) const
{
	return std::tie(axis, lo, size) < std::tie(other.axis, other.lo, other.size);
}

CartesianSubgrid::CartesianSubgrid(Index cells, Point spacing, Box refined, double epsilon, double permeability,
								   bool uniformFine, Index ratio)
	: CartesianSubgrid(cells, spacing, epsilon, permeability, {{refined, ratio}}, uniformFine)
{
}

CartesianSubgrid::CartesianSubgrid(Index cells, Point spacing, double epsilon, double permeability,
								   const std::vector<Region> &regions, bool uniformFine)
	: CartesianSubgrid(UniformLines(cells, spacing), epsilon, permeability, regions, uniformFine)
{
}

CartesianSubgrid::CartesianSubgrid(const Axes &lines, double epsilon, double permeability,
								   const std::vector<Region> &regions, bool uniformFine,
								   const MediumSampler &medium)
	: m_Ratio{1, 1, 1}, m_Epsilon(epsilon), m_Mu(permeability), m_UniformFine(uniformFine)
{
	Require(epsilon > 0 && std::isfinite(epsilon) && permeability > 0 && std::isfinite(permeability),
			"CartesianSubgrid: positive finite homogeneous epsilon and mu required");
	Require(!regions.empty() && regions.size() <= 64, "CartesianSubgrid: expected 1 through 64 regions");
	std::size_t count = 1;
	for (unsigned int n = 0; n < 3; ++n)
	{
		Require(lines[n].size() >= 5 && lines[n].size() <= 1025,
				"CartesianSubgrid: expected 4 through 1024 base intervals per axis");
		m_Cells[n] = unsigned(lines[n].size() - 1);
		for (double value : lines[n])
		{
			Require(std::isfinite(value), "CartesianSubgrid: finite base coordinates required");
		}
		for (unsigned int i = 0; i < m_Cells[n]; ++i)
		{
			const double width = lines[n][i + 1] - lines[n][i];
			Require(width > 0 && std::isfinite(width),
					"CartesianSubgrid: strictly increasing finite base intervals required");
			CheckSubdivision(lines[n][i], lines[n][i + 1], 1);
		}
		count *= m_Cells[n];
	}
	Require(count <= 1000000, "CartesianSubgrid: prototype allocation guard exceeds one million leaf cells");
	std::size_t leaves = count;
	for (std::size_t i = 0; i < regions.size(); ++i)
	{
		const Region &region = regions[i];
		ValidateRatio(region.ratio);
		std::size_t selected = 1;
		for (unsigned int n = 0; n < 3; ++n)
		{
			Require(region.box.start[n] >= 1 && region.box.start[n] < region.box.stop[n] &&
						region.box.stop[n] < m_Cells[n],
					"CartesianSubgrid: refinement box must have a base-cell margin on all six sides");
			m_Ratio[n] = std::max(m_Ratio[n], region.ratio[n]);
			for (unsigned int k = region.box.start[n]; k < region.box.stop[n]; ++k)
			{
				CheckSubdivision(lines[n][k], lines[n][k + 1], region.ratio[n]);
			}
			selected *= region.box.stop[n] - region.box.start[n];
		}
		for (std::size_t j = 0; j < i; ++j)
		{
			bool separated = false;
			for (unsigned int n = 0; n < 3; ++n)
			{
				separated = separated || region.box.stop[n] < regions[j].box.start[n] ||
							regions[j].box.stop[n] < region.box.start[n];
			}
			Require(separated, "CartesianSubgrid: regions require a full base-cell gap; touching, overlap and nesting "
							   "are unsupported");
		}
		const std::size_t children = std::size_t(region.ratio[0]) * region.ratio[1] * region.ratio[2];
		// Divide before multiplying; the accumulated budget covers all regions.
		Require(children - 1 <= (1000000 - leaves) / selected,
				"CartesianSubgrid: prototype allocation guard exceeds one million leaf cells");
		leaves += (children - 1) * selected;
	}
	const std::size_t fineChildren = std::size_t(m_Ratio[0]) * m_Ratio[1] * m_Ratio[2];
	Require(!uniformFine || count <= 1000000 / fineChildren,
			"CartesianSubgrid: prototype allocation guard exceeds one million leaf cells");
	if (uniformFine)
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int k = 0; k < m_Cells[n]; ++k)
			{
				CheckSubdivision(lines[n][k], lines[n][k + 1], m_Ratio[n]);
			}
		}
	}
	m_Regions = regions;
	m_Lines = lines;
	Build(medium);
}

void CartesianSubgrid::ValidateRatio(const Index &ratio)
{
	for (unsigned int value : ratio)
	{
		Require(value >= 1 && value <= 16 && (value & (value - 1)) == 0,
				"CartesianSubgrid: axis ratios must be powers of two from 1 to 16");
	}
}

std::size_t CartesianSubgrid::RootIndex(const Index &cell) const
{
	return (std::size_t(cell[0]) * m_Cells[1] + cell[1]) * m_Cells[2] + cell[2];
}

CartesianSubgrid::Edge CartesianSubgrid::LocalEdge(const Leaf &cell, unsigned int axis, unsigned int first,
												   unsigned int second) const
{
	Edge edge{int(axis), cell.lo, cell.size[axis]};
	edge.lo[(axis + 1) % 3] += first * cell.size[(axis + 1) % 3];
	edge.lo[(axis + 2) % 3] += second * cell.size[(axis + 2) % 3];
	return edge;
}

bool CartesianSubgrid::Boundary(const Edge &edge) const
{
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (n != unsigned(edge.axis) && (edge.lo[n] == 0 || edge.lo[n] == int(m_Ratio[n] * m_Cells[n])))
		{
			return true;
		}
	}
	return false;
}

bool CartesianSubgrid::Hanging(const Edge &edge) const
{
	return InterfaceBox(edge) != nullptr;
}

const CartesianSubgrid::Box *CartesianSubgrid::InterfaceBox(const Edge &edge) const
{
	if (m_UniformFine)
	{
		return nullptr;
	}
	bool coarse = edge.size == int(m_Ratio[edge.axis]);
	for (unsigned int n = 0; n < 3; ++n)
	{
		coarse = coarse && edge.lo[n] % int(m_Ratio[n]) == 0;
	}
	if (coarse)
	{
		return nullptr;
	}
	for (const Region &region : m_Regions)
	{
		bool inside = true, boundary = false;
		for (unsigned int n = 0; n < 3; ++n)
		{
			const int lo = int(m_Ratio[n] * region.box.start[n]);
			const int hi = int(m_Ratio[n] * region.box.stop[n]);
			inside = inside && edge.lo[n] >= lo && edge.lo[n] + (n == unsigned(edge.axis) ? edge.size : 0) <= hi;
			boundary = boundary || (n != unsigned(edge.axis) && (edge.lo[n] == lo || edge.lo[n] == hi));
		}
		if (inside && boundary)
		{
			return &region.box;
		}
	}
	return nullptr;
}

CartesianSubgrid::Trace CartesianSubgrid::Restrict(const Edge &edge) const
{
	Trace trace;
	if (Boundary(edge))
	{
		return trace;
	}
	const Box *box = InterfaceBox(edge);
	if (!box)
	{
		trace.weights.emplace_back(m_EdgeMap.at(edge), 1.);
		return trace;
	}
	// The coarse hexahedron's tangential edge shape is constant along the
	// edge direction and bilinear in the other two coordinates. On a face
	// this reduces to linear interpolation; at an edge/corner it reduces to
	// the single shared coarse edge. All six faces use this same rule.
	const unsigned int a = (edge.axis + 1) % 3;
	const unsigned int b = (edge.axis + 2) % 3;
	Index root;
	for (unsigned int n = 0; n < 3; ++n)
	{
		root[n] = std::min(unsigned(edge.lo[n] / int(m_Ratio[n])), m_Cells[n] - 1);
	}
	for (unsigned int n : {a, b})
	{
		if (edge.lo[n] == int(m_Ratio[n] * box->start[n]))
		{
			root[n] = box->start[n] - 1;
			break;
		}
		if (edge.lo[n] == int(m_Ratio[n] * box->stop[n]))
		{
			root[n] = box->stop[n];
			break;
		}
	}
	const double u = double(edge.lo[a]) / m_Ratio[a] - root[a];
	const double v = double(edge.lo[b]) / m_Ratio[b] - root[b];
	Require(u >= 0 && u <= 1 && v >= 0 && v <= 1, "CartesianSubgrid: invalid hanging-edge coordinates");
	Leaf coarse{{int(m_Ratio[0] * root[0]), int(m_Ratio[1] * root[1]), int(m_Ratio[2] * root[2])},
				{int(m_Ratio[0]), int(m_Ratio[1]), int(m_Ratio[2])}};
	for (unsigned int i = 0; i < 2; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			const double weight = (i ? u : 1 - u) * (j ? v : 1 - v);
			if (weight != 0)
			{
				const Edge master = LocalEdge(coarse, edge.axis, i, j);
				if (!Boundary(master))
				{
					trace.weights.emplace_back(m_EdgeMap.at(master), weight);
				}
			}
		}
	}
	return trace;
}

void CartesianSubgrid::Build(const MediumSampler &medium)
{
	const std::size_t roots = std::size_t(m_Cells[0]) * m_Cells[1] * m_Cells[2];
	m_RootLeaves.resize(roots);
	m_RootRatios.resize(roots, Index{1, 1, 1});
	for (unsigned int x = 0; x < m_Cells[0]; ++x)
	{
		for (unsigned int y = 0; y < m_Cells[1]; ++y)
		{
			for (unsigned int z = 0; z < m_Cells[2]; ++z)
			{
				const Index root{x, y, z};
				Index ratio{1, 1, 1};
				if (m_UniformFine)
				{
					ratio = m_Ratio;
				}
				else
				{
					for (const Region &region : m_Regions)
					{
						bool inside = true;
						for (unsigned int n = 0; n < 3; ++n)
						{
							inside = inside && root[n] >= region.box.start[n] && root[n] < region.box.stop[n];
						}
						if (inside)
						{
							ratio = region.ratio;
							break;
						}
					}
				}
				const std::size_t index = RootIndex(root);
				m_RootLeaves[index] = m_Leaves.size();
				m_RootRatios[index] = ratio;
				const std::array<int, 3> size{int(m_Ratio[0] / ratio[0]), int(m_Ratio[1] / ratio[1]),
											  int(m_Ratio[2] / ratio[2])};
				for (unsigned int i = 0; i < ratio[0]; ++i)
				{
					for (unsigned int j = 0; j < ratio[1]; ++j)
					{
						for (unsigned int k = 0; k < ratio[2]; ++k)
						{
							m_Leaves.push_back(
								{{int(m_Ratio[0] * x) + int(i) * size[0], int(m_Ratio[1] * y) + int(j) * size[1],
								  int(m_Ratio[2] * z) + int(k) * size[2]},
								 size});
						}
					}
				}
			}
		}
	}
	m_Media.reserve(m_Leaves.size());
	for (std::size_t index = 0; index < m_Leaves.size(); ++index)
	{
		const auto start = LeafStart(index), length = LeafLengths(index);
		Point stop;
		for (unsigned int n = 0; n < 3; ++n)
		{
			stop[n] = start[n] + length[n];
		}
		const Medium value = medium ? medium(start, stop) : Medium{m_Epsilon, m_Mu};
		Require(value.epsilon > 0 && std::isfinite(value.epsilon) && value.permeability > 0 &&
					std::isfinite(value.permeability),
				"CartesianSubgrid: positive finite leaf epsilon and permeability required");
		m_Media.push_back(value);
	}
	for (const Leaf &cell : m_Leaves)
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i < 2; ++i)
			{
				for (unsigned int j = 0; j < 2; ++j)
				{
					const Edge edge = LocalEdge(cell, n, i, j);
					if (!Boundary(edge) && !Hanging(edge) && m_EdgeMap.emplace(edge, m_Edges.size()).second)
					{
						m_Edges.push_back(edge);
					}
				}
			}
		}
	}
	std::vector<double> mass(m_Edges.size(), 0);
	std::map<std::pair<std::size_t, std::size_t>, double> off;
	m_RowStart.push_back(0);
	for (std::size_t cellIndex = 0; cellIndex < m_Leaves.size(); ++cellIndex)
	{
		const Leaf &cell = m_Leaves[cellIndex];
		const Medium &material = m_Media[cellIndex];
		Point length;
		for (unsigned int n = 0; n < 3; ++n)
		{
			length[n] = Length(cell, n);
		}
		const double volume = length[0] * length[1] * length[2];
		// Assemble P^T M_e P, not independent field interpolation. Keep the
		// off-diagonal interface mass; row-sum lumping here changes dispersion.
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int i = 0; i < 2; ++i)
			{
				for (unsigned int j = 0; j < 2; ++j)
				{
					const Trace trace = Restrict(LocalEdge(cell, n, i, j));
					for (const auto &a : trace.weights)
					{
						mass[a.first] += material.epsilon * volume / 4 * a.second * a.second;
						for (const auto &b : trace.weights)
						{
							if (a.first < b.first)
							{
								off[{a.first, b.first}] +=
									material.epsilon * volume / 4 * a.second * b.second;
							}
						}
					}
				}
			}
		}
		for (unsigned int n = 0; n < 3; ++n)
		{
			const unsigned int a = (n + 1) % 3;
			const unsigned int b = (n + 2) % 3;
			const double area = length[a] * length[b];
			for (unsigned int side = 0; side < 2; ++side)
			{
				std::map<std::size_t, double> row;
				auto edge = [&](unsigned int direction, unsigned int crossAxis, unsigned int high, double sign)
				{
					Edge local{int(direction), cell.lo, cell.size[direction]};
					local.lo[n] += side * cell.size[n];
					local.lo[crossAxis] += high * cell.size[crossAxis];
					for (const auto &term : Restrict(local).weights)
					{
						row[term.first] += sign * length[direction] * term.second;
					}
				};
				edge(b, a, 1, 1);
				edge(b, a, 0, -1);
				edge(a, b, 1, -1);
				edge(a, b, 0, 1);
				for (const auto &term : row)
				{
					if (term.second != 0)
					{
						m_Columns.push_back(term.first);
						m_Values.push_back(term.second);
					}
				}
				m_RowStart.push_back(m_Values.size());
				// One half-volume per face. Coincident faces are deliberately
				// separate quadrature histories, not a global fine H lattice.
				m_SqrtReluctance.push_back(std::sqrt(volume / (2 * material.permeability * area * area)));
			}
		}
	}
	for (double value : mass)
	{
		Require(value > 0 && std::isfinite(value), "CartesianSubgrid: non-positive electric mass");
		m_SqrtMass.push_back(std::sqrt(value));
	}
	const double inverseMassBound = FactorMass(off);
	std::vector<double> columns(m_Edges.size(), 0);
	double maxRow = 0;
	for (std::size_t row = 0; row + 1 < m_RowStart.size(); ++row)
	{
		double sum = 0;
		for (std::size_t i = m_RowStart[row]; i < m_RowStart[row + 1]; ++i)
		{
			m_Values[i] *= m_SqrtReluctance[row] / m_SqrtMass[m_Columns[i]];
			sum += std::abs(m_Values[i]);
			columns[m_Columns[i]] += std::abs(m_Values[i]);
		}
		maxRow = std::max(maxRow, sum);
	}
	const double bound = std::sqrt(maxRow * *std::max_element(columns.begin(), columns.end()) * inverseMassBound);
	Require(bound > 0 && std::isfinite(bound), "CartesianSubgrid: invalid coupled spectral bound");
	m_StableStep = 1.98 / bound;
	m_Electric.assign(m_Edges.size(), 0);
	m_Magnetic.assign(m_SqrtReluctance.size(), 0);
	m_Increment.resize(m_Electric.size());
	m_Curl.resize(m_Magnetic.size());
}

void CartesianSubgrid::MassBlock::Solve(std::vector<double> &values) const
{
	const std::size_t n = edges.size();
	for (std::size_t i = 0; i < n; ++i)
	{
		for (std::size_t j = 0; j < i; ++j)
		{
			values[edges[i]] -= lower[i * n + j] * values[edges[j]];
		}
		values[edges[i]] /= lower[i * n + i];
	}
	for (std::size_t i = n; i-- > 0;)
	{
		for (std::size_t j = i + 1; j < n; ++j)
		{
			values[edges[i]] -= lower[j * n + i] * values[edges[j]];
		}
		values[edges[i]] /= lower[i * n + i];
	}
}

double CartesianSubgrid::FactorMass(const std::map<std::pair<std::size_t, std::size_t>, double> &off)
{
	std::map<std::size_t, std::vector<std::size_t>> graph;
	for (const auto &term : off)
	{
		const std::size_t a = term.first.first, b = term.first.second;
		graph[a].push_back(b);
		graph[b].push_back(a);
		m_MassOffDiagonal.push_back({term.first, term.second / (m_SqrtMass[a] * m_SqrtMass[b])});
	}
	std::vector<bool> seen(m_Edges.size(), false);
	double inverseBound = 1;
	for (const auto &node : graph)
	{
		if (seen[node.first])
		{
			continue;
		}
		MassBlock block;
		block.edges.push_back(node.first);
		seen[node.first] = true;
		for (std::size_t i = 0; i < block.edges.size(); ++i)
		{
			for (std::size_t neighbour : graph[block.edges[i]])
			{
				if (!seen[neighbour])
				{
					seen[neighbour] = true;
					block.edges.push_back(neighbour);
				}
			}
		}
		const std::size_t n = block.edges.size();
		Require(n <= 1024, "CartesianSubgrid: prototype interface mass block exceeds 1024 edges");
		block.lower.assign(n * n, 0);
		for (std::size_t i = 0; i < n; ++i)
		{
			block.lower[i * n + i] = 1;
			for (std::size_t j = 0; j < i; ++j)
			{
				const auto pair = std::minmax(block.edges[i], block.edges[j]);
				const auto found = off.find({pair.first, pair.second});
				if (found != off.end())
				{
					block.lower[i * n + j] = found->second / (m_SqrtMass[pair.first] * m_SqrtMass[pair.second]);
				}
			}
		}
		for (std::size_t i = 0; i < n; ++i)
		{
			for (std::size_t j = 0; j <= i; ++j)
			{
				double value = block.lower[i * n + j];
				for (std::size_t k = 0; k < j; ++k)
				{
					value -= block.lower[i * n + k] * block.lower[j * n + k];
				}
				if (i == j)
				{
					Require(value > 0 && std::isfinite(value), "CartesianSubgrid: interface mass is not SPD");
					block.lower[i * n + j] = std::sqrt(value);
				}
				else
				{
					block.lower[i * n + j] = value / block.lower[j * n + j];
				}
			}
		}
		// A conservative ||S^-1||_2 <= ||S^-1||_infinity bound, computed once.
		std::vector<double> rowSum(n, 0), column(m_Edges.size(), 0);
		for (std::size_t j = 0; j < n; ++j)
		{
			for (std::size_t index : block.edges)
			{
				column[index] = 0;
			}
			column[block.edges[j]] = 1;
			block.Solve(column);
			for (std::size_t i = 0; i < n; ++i)
			{
				rowSum[i] += std::abs(column[block.edges[i]]);
			}
		}
		inverseBound = std::max(inverseBound, *std::max_element(rowSum.begin(), rowSum.end()));
		m_Blocks.push_back(std::move(block));
	}
	return inverseBound * (1 + 1e-12);
}

void CartesianSubgrid::Curl(const std::vector<double> &electric, std::vector<double> &magnetic) const
{
	for (std::size_t row = 0; row < magnetic.size(); ++row)
	{
		double value = 0;
		for (std::size_t i = m_RowStart[row]; i < m_RowStart[row + 1]; ++i)
		{
			value += m_Values[i] * electric[m_Columns[i]];
		}
		magnetic[row] = value;
	}
}

void CartesianSubgrid::CurlTranspose(const std::vector<double> &magnetic, std::vector<double> &electric) const
{
	std::fill(electric.begin(), electric.end(), 0);
	for (std::size_t row = 0; row < magnetic.size(); ++row)
	{
		for (std::size_t i = m_RowStart[row]; i < m_RowStart[row + 1]; ++i)
		{
			electric[m_Columns[i]] += m_Values[i] * magnetic[row];
		}
	}
}

void CartesianSubgrid::SolveMass(std::vector<double> &values) const
{
	for (const MassBlock &block : m_Blocks)
	{
		block.Solve(values);
	}
}

void CartesianSubgrid::CheckStep(double timestep) const
{
	Require(timestep > 0 && std::isfinite(timestep) && timestep <= m_StableStep,
			"CartesianSubgrid: timestep exceeds the coupled stability bound");
}

void CartesianSubgrid::SetSources(const std::vector<std::pair<std::size_t, double>> &sources)
{
	m_Sources.clear();
	for (const auto &source : sources)
	{
		Require(source.first < m_Electric.size() && std::isfinite(source.second),
				"CartesianSubgrid: invalid source edge");
		m_Sources.emplace_back(source.first, source.second * m_SqrtMass[source.first]);
	}
}

void CartesianSubgrid::Step(double timestep, double sourceAmplitude)
{
	CheckStep(timestep);
	Require(std::isfinite(sourceAmplitude), "CartesianSubgrid: nonfinite drive");
	CurlTranspose(m_Magnetic, m_Increment);
	SolveMass(m_Increment);
	for (std::size_t i = 0; i < m_Electric.size(); ++i)
	{
		m_Electric[i] += timestep * m_Increment[i];
	}
	for (const auto &source : m_Sources)
	{
		m_Electric[source.first] += sourceAmplitude * source.second;
	}
	Curl(m_Electric, m_Curl);
	for (std::size_t i = 0; i < m_Magnetic.size(); ++i)
	{
		m_Magnetic[i] -= timestep * m_Curl[i];
	}
}

double CartesianSubgrid::Coordinate(unsigned int axis, double lattice) const
{
	const unsigned int root = std::min(unsigned(lattice / m_Ratio[axis]), m_Cells[axis] - 1);
	const double fraction = (lattice - double(root) * m_Ratio[axis]) / m_Ratio[axis];
	if (fraction == 0)
	{
		return m_Lines[axis][root];
	}
	if (fraction == 1)
	{
		return m_Lines[axis][root + 1];
	}
	return m_Lines[axis][root] + fraction * (m_Lines[axis][root + 1] - m_Lines[axis][root]);
}

double CartesianSubgrid::Length(const Leaf &leaf, unsigned int axis) const
{
	const unsigned int root = unsigned(leaf.lo[axis]) / m_Ratio[axis];
	return (m_Lines[axis][root + 1] - m_Lines[axis][root]) * leaf.size[axis] / m_Ratio[axis];
}

CartesianSubgrid::Point CartesianSubgrid::LeafStart(std::size_t leaf) const
{
	const auto &cell = m_Leaves.at(leaf);
	return {Coordinate(0, cell.lo[0]), Coordinate(1, cell.lo[1]), Coordinate(2, cell.lo[2])};
}

CartesianSubgrid::Point CartesianSubgrid::LeafLengths(std::size_t leaf) const
{
	const auto &cell = m_Leaves.at(leaf);
	return {Length(cell, 0), Length(cell, 1), Length(cell, 2)};
}

CartesianSubgrid::Point CartesianSubgrid::EdgePosition(std::size_t edge) const
{
	const Edge &key = m_Edges.at(edge);
	Point position;
	for (unsigned int n = 0; n < 3; ++n)
	{
		position[n] = Coordinate(n, key.lo[n] + (n == unsigned(key.axis) ? key.size / 2. : 0));
	}
	return position;
}

unsigned int CartesianSubgrid::EdgeDirection(std::size_t edge) const
{
	return m_Edges.at(edge).axis;
}

void CartesianSubgrid::SetElectric(const Field &field)
{
	for (std::size_t edge = 0; edge < m_Electric.size(); ++edge)
	{
		const double value = field(EdgeDirection(edge), EdgePosition(edge));
		Require(std::isfinite(value), "CartesianSubgrid: nonfinite initial electric field");
		m_Electric[edge] = value * m_SqrtMass[edge];
	}
}

void CartesianSubgrid::SetMagnetic(const Field &field)
{
	for (std::size_t cell = 0; cell < m_Leaves.size(); ++cell)
	{
		const Leaf &leaf = m_Leaves[cell];
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int side = 0; side < 2; ++side)
			{
				Point position;
				for (unsigned int k = 0; k < 3; ++k)
				{
					position[k] = Coordinate(k, leaf.lo[k] + leaf.size[k] * (n == k ? side : .5));
				}
				const double area = Length(leaf, (n + 1) % 3) * Length(leaf, (n + 2) % 3);
				const std::size_t index = 6 * cell + 2 * n + side;
				const double value = field(n, position);
				Require(std::isfinite(value), "CartesianSubgrid: nonfinite initial magnetic field");
				m_Magnetic[index] = value * m_Media[cell].permeability * area * m_SqrtReluctance[index];
			}
		}
	}
}

std::size_t CartesianSubgrid::FindLeaf(const Point &position) const
{
	Index root;
	for (unsigned int n = 0; n < 3; ++n)
	{
		Require(std::isfinite(position[n]) && position[n] >= m_Lines[n].front() && position[n] <= m_Lines[n].back(),
				"CartesianSubgrid: sample outside domain");
		const auto next = std::upper_bound(m_Lines[n].begin(), m_Lines[n].end(), position[n]);
		root[n] = std::min(unsigned(next - m_Lines[n].begin() - 1), m_Cells[n] - 1);
		const double width = m_Lines[n][root[n] + 1] - m_Lines[n][root[n]];
		const double boundary = m_Lines[n][root[n] + 1];
		const double tolerance = std::min(width / m_Ratio[n] * 1e-8,
										  32 * std::numeric_limits<double>::epsilon() *
											  std::max({std::abs(position[n]), std::abs(boundary), width}));
		// Material-normal fields can jump. At a cell face select the positive
		// side consistently despite SI/drawing-unit coordinate roundoff.
		// Only lookup is disambiguated; geometry and sample positions stay fixed.
		if (root[n] + 1 < m_Cells[n] && std::abs(position[n] - boundary) <= tolerance)
		{
			++root[n];
		}
	}
	const std::size_t index = RootIndex(root);
	const std::size_t first = m_RootLeaves[index];
	const Index &ratio = m_RootRatios[index];
	Index child;
	for (unsigned int n = 0; n < 3; ++n)
	{
		const double start = m_Lines[n][root[n]];
		const double width = m_Lines[n][root[n] + 1] - start;
		double local = ratio[n] * (position[n] - start) / width;
		const double nearest = std::clamp(std::round(local), 0., double(ratio[n]));
		const double boundary = start + nearest * width / ratio[n];
		const double tolerance = std::min(width / ratio[n] * 1e-8,
										  32 * std::numeric_limits<double>::epsilon() *
											  std::max({std::abs(position[n]), std::abs(boundary), width}));
		if (std::abs(position[n] - boundary) <= tolerance)
		{
			local = nearest;
		}
		child[n] = std::min(unsigned(std::max(0., local)), ratio[n] - 1);
	}
	return first + (std::size_t(child[0]) * ratio[1] + child[1]) * ratio[2] + child[2];
}

double CartesianSubgrid::TraceValue(const Trace &trace, const std::vector<double> &amplitude) const
{
	double value = 0;
	for (const auto &term : trace.weights)
	{
		value += term.second * amplitude[term.first] / m_SqrtMass[term.first];
	}
	return value;
}

double CartesianSubgrid::Electric(unsigned int n, const Point &position) const
{
	Require(n < 3, "CartesianSubgrid: invalid component");
	const Leaf &cell = m_Leaves[FindLeaf(position)];
	const unsigned int a = (n + 1) % 3, b = (n + 2) % 3;
	const double u = (position[a] - Coordinate(a, cell.lo[a])) / Length(cell, a);
	const double v = (position[b] - Coordinate(b, cell.lo[b])) / Length(cell, b);
	double value = 0;
	for (unsigned int i = 0; i < 2; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			value += (i ? u : 1 - u) * (j ? v : 1 - v) * TraceValue(Restrict(LocalEdge(cell, n, i, j)), m_Electric);
		}
	}
	return value;
}

double CartesianSubgrid::Magnetic(unsigned int n, const Point &position) const
{
	Require(n < 3, "CartesianSubgrid: invalid component");
	const std::size_t index = FindLeaf(position);
	const Leaf &cell = m_Leaves[index];
	const double t = (position[n] - Coordinate(n, cell.lo[n])) / Length(cell, n);
	const double area = Length(cell, (n + 1) % 3) * Length(cell, (n + 2) % 3);
	const std::size_t first = 6 * index + 2 * n;
	return ((1 - t) * m_Magnetic[first] / m_SqrtReluctance[first] +
			t * m_Magnetic[first + 1] / m_SqrtReluctance[first + 1]) /
		   (area * m_Media[index].permeability);
}

void CartesianSubgrid::AddElectric(std::size_t edge, double value)
{
	Require(std::isfinite(value), "CartesianSubgrid: nonfinite source");
	m_Electric.at(edge) += value * m_SqrtMass.at(edge);
}

double CartesianSubgrid::ElectricEnergy() const
{
	double result = .5 * std::inner_product(m_Electric.begin(), m_Electric.end(), m_Electric.begin(), 0.);
	for (const auto &term : m_MassOffDiagonal)
	{
		result += term.second * m_Electric[term.first.first] * m_Electric[term.first.second];
	}
	return result;
}

double CartesianSubgrid::MagneticEnergy() const
{
	return .5 * std::inner_product(m_Magnetic.begin(), m_Magnetic.end(), m_Magnetic.begin(), 0.);
}

double CartesianSubgrid::Energy(double timestep) const
{
	double value = ElectricEnergy() + MagneticEnergy();
	if (timestep != 0)
	{
		CheckStep(timestep);
		std::vector<double> curl(m_Magnetic.size());
		Curl(m_Electric, curl);
		value += .5 * timestep * std::inner_product(curl.begin(), curl.end(), m_Magnetic.begin(), 0.);
	}
	return value;
}

double CartesianSubgrid::CurlDivergenceResidual(const Field &field) const
{
	std::vector<double> electric(m_Edges.size()), curl(m_Magnetic.size());
	for (std::size_t i = 0; i < electric.size(); ++i)
	{
		electric[i] = field(EdgeDirection(i), EdgePosition(i)) * m_SqrtMass[i];
	}
	Curl(electric, curl);
	double maximum = 0;
	for (std::size_t cell = 0; cell < m_Leaves.size(); ++cell)
	{
		double sum = 0, scale = 0;
		for (unsigned int n = 0; n < 3; ++n)
		{
			const std::size_t i = 6 * cell + 2 * n;
			sum += curl[i + 1] / m_SqrtReluctance[i + 1] - curl[i] / m_SqrtReluctance[i];
			scale += std::abs(curl[i + 1] / m_SqrtReluctance[i + 1]) + std::abs(curl[i] / m_SqrtReluctance[i]);
		}
		maximum = std::max(maximum, std::abs(sum) / std::max(scale, 1e-300));
	}
	return maximum;
}

std::size_t CartesianSubgrid::StorageBytes() const
{
	std::size_t bytes =
		sizeof(*this) + m_Leaves.capacity() * sizeof(Leaf) + m_RootLeaves.capacity() * sizeof(std::size_t);
	for (const auto &axis : m_Lines)
	{
		bytes += axis.capacity() * sizeof(double);
	}
	bytes += m_Edges.capacity() * sizeof(Edge) + m_RootRatios.capacity() * sizeof(Index) +
			 m_Regions.capacity() * sizeof(Region) + m_Media.capacity() * sizeof(Medium);
	// std::map node overhead is implementation-specific; this is an approximate
	// accounting estimate, not a measurement of allocator or process RSS.
	bytes += m_EdgeMap.size() * (sizeof(Edge) + sizeof(std::size_t) + 4 * sizeof(void *));
	for (const auto *values :
		 {&m_SqrtMass, &m_SqrtReluctance, &m_Values, &m_Electric, &m_Magnetic, &m_Increment, &m_Curl})
	{
		bytes += values->capacity() * sizeof(double);
	}
	bytes += (m_RowStart.capacity() + m_Columns.capacity()) * sizeof(std::size_t);
	bytes += m_MassOffDiagonal.capacity() * sizeof(decltype(m_MassOffDiagonal)::value_type);
	bytes += m_Sources.capacity() * sizeof(decltype(m_Sources)::value_type);
	for (const MassBlock &block : m_Blocks)
	{
		bytes += sizeof(block) + block.edges.capacity() * sizeof(std::size_t) + block.lower.capacity() * sizeof(double);
	}
	return bytes;
}

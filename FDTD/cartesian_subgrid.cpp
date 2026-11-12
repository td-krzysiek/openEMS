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
} // namespace

bool CartesianSubgrid::Edge::operator<(const Edge &other) const
{
	return std::tie(axis, lo, size) < std::tie(other.axis, other.lo, other.size);
}

CartesianSubgrid::CartesianSubgrid(Index cells, Point spacing, Box refined, double epsilon, double permeability,
								   bool uniformFine)
	: m_Cells(cells), m_Spacing(spacing), m_Box(refined), m_Epsilon(epsilon), m_Mu(permeability),
	  m_UniformFine(uniformFine)
{
	Require(epsilon > 0 && std::isfinite(epsilon) && permeability > 0 && std::isfinite(permeability),
			"CartesianSubgrid: positive finite homogeneous epsilon and mu required");
	std::size_t count = 1;
	for (unsigned int n = 0; n < 3; ++n)
	{
		Require(cells[n] >= 4 && cells[n] <= 1024 && spacing[n] > 0 && std::isfinite(spacing[n]),
				"CartesianSubgrid: invalid base dimensions or spacing");
		Require(refined.start[n] >= 1 && refined.start[n] < refined.stop[n] && refined.stop[n] < cells[n],
				"CartesianSubgrid: refinement box must have a base-cell margin on all six sides");
		count *= cells[n];
	}
	Require(count <= 20000000, "CartesianSubgrid: prototype base-cell limit exceeded");
	const std::size_t refinedCells = std::size_t(refined.stop[0] - refined.start[0]) *
									 (refined.stop[1] - refined.start[1]) * (refined.stop[2] - refined.start[2]);
	Require((uniformFine ? 8 * count : count + 7 * refinedCells) <= 1000000,
			"CartesianSubgrid: prototype allocation guard exceeds one million leaf cells");
	Build();
}

std::size_t CartesianSubgrid::RootIndex(const Index &cell) const
{
	return (std::size_t(cell[0]) * m_Cells[1] + cell[1]) * m_Cells[2] + cell[2];
}

CartesianSubgrid::Edge CartesianSubgrid::LocalEdge(const Leaf &cell, unsigned int axis, unsigned int first,
												   unsigned int second) const
{
	Edge edge{int(axis), cell.lo, cell.size};
	edge.lo[(axis + 1) % 3] += first * cell.size;
	edge.lo[(axis + 2) % 3] += second * cell.size;
	return edge;
}

bool CartesianSubgrid::Boundary(const Edge &edge) const
{
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (n != unsigned(edge.axis) && (edge.lo[n] == 0 || edge.lo[n] == int(2 * m_Cells[n])))
		{
			return true;
		}
	}
	return false;
}

bool CartesianSubgrid::Hanging(const Edge &edge) const
{
	if (m_UniformFine || edge.size != 1)
	{
		return false;
	}
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (n != unsigned(edge.axis) && (edge.lo[n] == int(2 * m_Box.start[n]) || edge.lo[n] == int(2 * m_Box.stop[n])))
		{
			return true;
		}
	}
	return false;
}

CartesianSubgrid::Trace CartesianSubgrid::Restrict(const Edge &edge) const
{
	Trace trace;
	if (Boundary(edge))
	{
		return trace;
	}
	if (!Hanging(edge))
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
		root[n] = std::min(unsigned(edge.lo[n] / 2), m_Cells[n] - 1);
	}
	for (unsigned int n : {a, b})
	{
		if (edge.lo[n] == int(2 * m_Box.start[n]))
		{
			root[n] = m_Box.start[n] - 1;
			break;
		}
		if (edge.lo[n] == int(2 * m_Box.stop[n]))
		{
			root[n] = m_Box.stop[n];
			break;
		}
	}
	const double u = edge.lo[a] / 2. - root[a];
	const double v = edge.lo[b] / 2. - root[b];
	Require(u >= 0 && u <= 1 && v >= 0 && v <= 1, "CartesianSubgrid: invalid hanging-edge coordinates");
	Leaf coarse{{int(2 * root[0]), int(2 * root[1]), int(2 * root[2])}, 2};
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

void CartesianSubgrid::Build()
{
	m_RootLeaves.resize(std::size_t(m_Cells[0]) * m_Cells[1] * m_Cells[2]);
	for (unsigned int x = 0; x < m_Cells[0]; ++x)
	{
		for (unsigned int y = 0; y < m_Cells[1]; ++y)
		{
			for (unsigned int z = 0; z < m_Cells[2]; ++z)
			{
				const Index root{x, y, z};
				bool refined = m_UniformFine;
				bool inside = true;
				for (unsigned int n = 0; n < 3; ++n)
				{
					inside = inside && root[n] >= m_Box.start[n] && root[n] < m_Box.stop[n];
				}
				refined = refined || inside;
				m_RootLeaves[RootIndex(root)] = m_Leaves.size();
				for (int i = 0; i < (refined ? 2 : 1); ++i)
				{
					for (int j = 0; j < (refined ? 2 : 1); ++j)
					{
						for (int k = 0; k < (refined ? 2 : 1); ++k)
						{
							m_Leaves.push_back({{int(2 * x) + i, int(2 * y) + j, int(2 * z) + k}, refined ? 1 : 2});
						}
					}
				}
			}
		}
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
	for (const Leaf &cell : m_Leaves)
	{
		Point length;
		for (unsigned int n = 0; n < 3; ++n)
		{
			length[n] = m_Spacing[n] * cell.size / 2.;
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
						mass[a.first] += m_Epsilon * volume / 4 * a.second * a.second;
						for (const auto &b : trace.weights)
						{
							if (a.first < b.first)
							{
								off[{a.first, b.first}] += m_Epsilon * volume / 4 * a.second * b.second;
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
					Edge local{int(direction), cell.lo, cell.size};
					local.lo[n] += side * cell.size;
					local.lo[crossAxis] += high * cell.size;
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
				m_SqrtReluctance.push_back(std::sqrt(volume / (2 * m_Mu * area * area)));
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

CartesianSubgrid::Point CartesianSubgrid::EdgePosition(std::size_t edge) const
{
	const Edge &key = m_Edges.at(edge);
	Point position;
	for (unsigned int n = 0; n < 3; ++n)
	{
		position[n] = m_Spacing[n] * (key.lo[n] / 2. + (n == unsigned(key.axis) ? key.size / 4. : 0));
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
					position[k] = m_Spacing[k] * (leaf.lo[k] / 2. + leaf.size / 2. * (n == k ? side : .5));
				}
				const double area = m_Spacing[(n + 1) % 3] * m_Spacing[(n + 2) % 3] * leaf.size * leaf.size / 4.;
				const std::size_t index = 6 * cell + 2 * n + side;
				const double value = field(n, position);
				Require(std::isfinite(value), "CartesianSubgrid: nonfinite initial magnetic field");
				m_Magnetic[index] = value * m_Mu * area * m_SqrtReluctance[index];
			}
		}
	}
}

std::size_t CartesianSubgrid::FindLeaf(const Point &position) const
{
	Index root;
	for (unsigned int n = 0; n < 3; ++n)
	{
		Require(std::isfinite(position[n]) && position[n] >= 0 && position[n] <= m_Spacing[n] * m_Cells[n],
				"CartesianSubgrid: sample outside domain");
		root[n] = std::min(unsigned(position[n] / m_Spacing[n]), m_Cells[n] - 1);
	}
	const std::size_t first = m_RootLeaves[RootIndex(root)];
	if (m_Leaves[first].size == 2)
	{
		return first;
	}
	Index child;
	for (unsigned int n = 0; n < 3; ++n)
	{
		child[n] = std::min(unsigned(2 * (position[n] / m_Spacing[n] - root[n])), 1u);
	}
	return first + 4 * child[0] + 2 * child[1] + child[2];
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
	const double u = (position[a] / m_Spacing[a] * 2 - cell.lo[a]) / cell.size;
	const double v = (position[b] / m_Spacing[b] * 2 - cell.lo[b]) / cell.size;
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
	const double t = (position[n] / m_Spacing[n] * 2 - cell.lo[n]) / cell.size;
	const double area = m_Spacing[(n + 1) % 3] * m_Spacing[(n + 2) % 3] * cell.size * cell.size / 4.;
	const std::size_t first = 6 * index + 2 * n;
	return ((1 - t) * m_Magnetic[first] / m_SqrtReluctance[first] +
			t * m_Magnetic[first + 1] / m_SqrtReluctance[first + 1]) /
		   (area * m_Mu);
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
	bytes += m_Edges.capacity() * sizeof(Edge);
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

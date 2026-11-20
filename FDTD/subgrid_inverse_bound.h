/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_INVERSE_BOUND_H
#define SUBGRID_INVERSE_BOUND_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace CartesianSubgridDetail
{
//! Exact inverse infinity norm for an SPD matrix with nonnegative off-diagonals.
/*!
 * The caller supplies its undirected nonzero graph and an existing factor solve.
 * If the graph is bipartite, a diagonal sign matrix D makes D S D Stieltjes:
 * its inverse is entrywise nonnegative. Hence abs(S^-1) 1 = D S^-1 D 1,
 * so a single signed RHS gives the same row sums as all inverse columns.
 * An odd cycle returns false, leaving the caller's general fallback in charge.
 * No matrix entry, factorization or iterative tolerance is changed.
 */
template <class Solve>
bool BipartiteInverseNorm(const std::vector<std::vector<std::size_t>> &graph, Solve solve, double &bound)
{
	std::vector<int> sign(graph.size(), 0);
	std::vector<std::size_t> queue;
	queue.reserve(graph.size());
	for (std::size_t root = 0; root < graph.size(); ++root)
	{
		if (sign[root])
		{
			continue;
		}
		sign[root] = 1;
		queue.clear();
		queue.push_back(root);
		for (std::size_t first = 0; first < queue.size(); ++first)
		{
			const std::size_t node = queue[first];
			for (std::size_t neighbour : graph[node])
			{
				if (neighbour >= graph.size())
				{
					throw std::invalid_argument("CartesianSubgrid: invalid inverse-bound graph");
				}
				if (!sign[neighbour])
				{
					sign[neighbour] = -sign[node];
					queue.push_back(neighbour);
				}
				else if (sign[neighbour] == sign[node])
				{
					return false;
				}
			}
		}
	}
	std::vector<double> rowSum(sign.begin(), sign.end());
	solve(rowSum);
	bound = 0;
	for (std::size_t i = 0; i < rowSum.size(); ++i)
	{
		const double value = sign[i] * rowSum[i];
		if (!(value > 0) || !std::isfinite(value))
		{
			return false;
		}
		bound = std::max(bound, value);
	}
	return true;
}
} // namespace CartesianSubgridDetail
#endif

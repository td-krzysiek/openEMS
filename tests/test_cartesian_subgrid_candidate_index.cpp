/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/subgrid_candidate_index.h"
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>

namespace
{
using Index = CartesianSubgridCandidateIndex<std::size_t>;
using Point = Index::Point;
using Bounds = Index::Bounds;

void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

bool Contains(const Bounds &box, const Point &point, double tolerance)
{
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (point[n] < box[2 * n] - tolerance || point[n] > box[2 * n + 1] + tolerance)
		{
			return false;
		}
	}
	return true;
}

bool Overlaps(const Bounds &box, const Point &lo, const Point &hi)
{
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (box[2 * n] > hi[n] || box[2 * n + 1] < lo[n])
		{
			return false;
		}
	}
	return true;
}
}

int main()
{
	try
	{
		const Bounds domain{-8, 8, -4, 4, -2, 2};
		std::vector<Bounds> bounds;
		for (unsigned int i = 0; i < 80; ++i)
		{
			const double x = -7.8 + 1.55 * (i % 10), y = -3.7 + .95 * (i / 10);
			bounds.push_back({x, x + .8, y, y + .7, -.1, .1});
		}
		// Nested/overlapping volumes, a plane, an outside object and a large
		// global object retain the caller's existing priority/insertion order.
		bounds.insert(bounds.begin() + 3, {-1, 1, -1, 1, -1, 1});
		bounds.insert(bounds.begin() + 4, {-.5, .5, -.5, .5, -.5, .5});
		bounds.push_back({-8, 8, -4, 4, 0, 0});
		bounds.push_back({50, 51, 50, 51, 50, 51});
		bounds.push_back(domain);
		std::vector<std::size_t> ids(bounds.size());
		std::iota(ids.begin(), ids.end(), 0);
		Index index(ids, bounds, domain, .0001);
		Require(index.Indexed(), "ordinary control did not build an index");
		std::mt19937 random(170070);
		std::uniform_real_distribution<double> uniform(0, 1);
		std::size_t checks = 0, candidates = 0;
		for (unsigned int trial = 0; trial < 100000; ++trial)
		{
			Point point, lo, hi;
			for (unsigned int n = 0; n < 3; ++n)
			{
				point[n] = domain[2 * n] + uniform(random) * (domain[2 * n + 1] - domain[2 * n]);
				lo[n] = point[n] - .003;
				hi[n] = point[n] + .003;
			}
			const auto &nearby = index.At(point.data());
			candidates += nearby.size();
			Require(std::is_sorted(nearby.begin(), nearby.end()), "candidate priority order changed");
			std::size_t original = ids.size(), selected = ids.size();
			for (std::size_t id : ids)
			{
				if (Contains(bounds[id], point, .0001))
				{
					original = id;
					break;
				}
			}
			for (std::size_t id : nearby)
			{
				if (Contains(bounds[id], point, .0001))
				{
					selected = id;
					break;
				}
			}
			Require(original == selected, "point broad phase changed the first matching object");
			const auto &overlap = index.Overlapping(lo.data(), hi.data());
			for (std::size_t id : ids)
			{
				if (Overlaps(bounds[id], lo, hi))
				{
					Require(std::find(overlap.begin(), overlap.end(), id) != overlap.end(),
						"leaf broad phase removed an intersecting bound");
				}
			}
			++checks;
		}
		Require(candidates < 100000 * ids.size() / 4, "localized fixture did not reduce candidates");
		for (const Bounds &box : bounds)
		{
			for (unsigned int corner = 0; corner < 8; ++corner)
			{
				Point point;
				for (unsigned int n = 0; n < 3; ++n)
				{
					point[n] = box[2 * n + ((corner >> n) & 1)];
				}
				const auto &nearby = index.At(point.data());
				for (std::size_t id : ids)
				{
					if (Contains(bounds[id], point, .0001))
					{
						Require(std::find(nearby.begin(), nearby.end(), id) != nearby.end(),
							"boundary or outside query lost a candidate");
					}
				}
				++checks;
			}
		}
		const Point nonfinite{std::numeric_limits<double>::quiet_NaN(), 0, 0};
		Require(index.At(nonfinite.data()) == ids && index.At(nullptr) == ids,
			"invalid query must retain original full list");
		const Point lo{-7, -3, -1}, hi{7, 3, 1};
		Require(index.Overlapping(lo.data(), hi.data()) == ids, "multi-bin leaf must reuse full list");
		Index capped(ids, bounds, domain, 0, 16, 1);
		Require(!capped.Indexed() && capped.Entries() == 0 && capped.At(lo.data()) == ids,
			"entry-budget fallback retained partial buckets");
		auto invalidBounds = bounds;
		invalidBounds[0][0] = std::numeric_limits<double>::infinity();
		Index invalid(ids, invalidBounds, domain, 0);
		Require(!invalid.Indexed() && invalid.At(lo.data()) == ids, "invalid bounds must retain full list");
		Bounds extreme = domain;
		extreme[0] = -std::numeric_limits<double>::max();
		extreme[1] = std::numeric_limits<double>::max();
		Index overflow(ids, bounds, extreme, 0);
		Require(!overflow.Indexed(), "overflowing domain extent must not enter bucket arithmetic");
		auto reverse = ids;
		std::reverse(reverse.begin(), reverse.end());
		auto reverseBounds = bounds;
		std::reverse(reverseBounds.begin(), reverseBounds.end());
		Index reversed(reverse, reverseBounds, domain, 0);
		const Point centre{0, 0, 0};
		const auto &order = reversed.At(centre.data());
		Require(std::is_sorted(order.rbegin(), order.rend()), "input order was replaced by object ID order");
		std::cout << "CANDIDATE_INDEX exact_queries=" << checks
			<< " bounded_fallback=1 priority_order=1 entries=" << index.Entries()
			<< " average_candidates=" << double(candidates) / 100000 << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

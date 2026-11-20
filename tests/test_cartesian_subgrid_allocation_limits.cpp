/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/subgrid_allocation_limits.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace D = CartesianSubgridDetail;
using Grid = CartesianSubgrid;

static void Require(bool ok, const char *message)
{
	if (!ok)
	{
		throw std::runtime_error(message);
	}
}

template <class Action> static void Reject(Action action)
{
	bool rejected = false;
	try
	{
		action();
	}
	catch (const std::invalid_argument &)
	{
		rejected = true;
	}
	Require(rejected, "invalid cell-count policy was accepted");
}

static Grid Fixture(std::size_t limit = D::DefaultLeafLimit)
{
	Grid::Axes axes;
	for (auto &axis : axes)
	{
		for (unsigned int k = 0; k <= 6; ++k)
		{
			axis.push_back(.001 * k);
		}
	}
	return Grid(axes, 8.85418781762e-12, 1.2566370614359173e-6,
		std::vector<Grid::Region>{{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}},
		false, {}, {}, {}, {}, Grid::CurlLayout::Patterns, {}, limit);
}

int main()
{
	try
	{
		Require(D::DefaultLeafLimit == 1000000 && D::MaximumIndexLeafLimit == 357913941,
			"default or device representability boundary changed");
		D::ValidateLeafLimit(1);
		D::ValidateLeafLimit(D::MaximumIndexLeafLimit);
		Reject([] { D::ValidateLeafLimit(0); });
		Reject([] { D::ValidateLeafLimit(D::MaximumIndexLeafLimit + 1); });
		Require(D::LeafProduct(123, 7, 861) == 861, "exact product boundary rejected");
		Require(D::RefinedLeafCount(216, 8, 8, 272) == 272, "exact refined budget rejected");
		Reject([] { D::LeafProduct(123, 7, 860); });
		Reject([] { D::LeafProduct(1, 0, 100); });
		Reject([] { D::LeafProduct(std::numeric_limits<std::size_t>::max(), 2, D::MaximumIndexLeafLimit); });
		Reject([] { D::LeafProduct(1, std::numeric_limits<std::size_t>::max(), D::MaximumIndexLeafLimit); });
		Reject([] { D::RefinedLeafCount(216, 8, 8, 271); });
		Reject([] { D::RefinedLeafCount(273, 8, 8, 272); });
		Reject([] { D::RefinedLeafCount(0, 1, std::numeric_limits<std::size_t>::max(), D::MaximumIndexLeafLimit); });
		Reject([] { D::RefinedLeafCount(0, std::numeric_limits<std::size_t>::max(), 2, D::MaximumIndexLeafLimit); });
		Reject([] { D::RefinedLeafCount(1, 0, 1, 10); });
		Reject([] { D::RefinedLeafCount(1, 1, 0, 10); });
		Reject([] { Fixture(271); });
		auto original = Fixture(), exact = Fixture(272), larger = Fixture(4000000);
		Require(original.CellCount() == 272 && exact.CellCount() == 272 && larger.CellCount() == 272,
			"cell-count option changed physical topology");
		for (auto *grid : {&original, &exact, &larger})
		{
			grid->SetElectric([](unsigned int n, const Grid::Point &p)
			{
				return std::sin(517 * p[0] + 943 * p[1] - 617 * p[2] + n);
			});
		}
		const double dt = .7 * original.StableTimestep();
		Require(exact.StableTimestep() == original.StableTimestep() &&
			larger.StableTimestep() == original.StableTimestep(), "cell policy changed the coupled bound");
		for (unsigned int k = 0; k < 256; ++k)
		{
			original.Step(dt);
			exact.Step(dt);
			larger.Step(dt);
		}
		const auto a = original.DeviceStateData(), b = exact.DeviceStateData(), c = larger.DeviceStateData();
		Require(a.electric == b.electric && a.electric == c.electric &&
			a.magnetic == b.magnetic && a.magnetic == c.magnetic,
			"cell-count option changed deterministic fields");
		std::cout << "LEAF_BUDGET actual=272 exact=272 default=1000000 opt_in=4000000 exact_fields=1 overflow_guards=1\n";
		std::cout << "PASS: explicit leaf budgets, preallocation rejection and unchanged fields\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

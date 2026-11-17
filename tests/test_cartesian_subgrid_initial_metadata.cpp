/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_child_knots_fixture.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>

using Plan = CartesianSubgridDomainPlan;

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

template <class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const std::logic_error &)
	{
		return;
	}
	throw std::runtime_error("invalid initial metadata accepted");
}

template <class Value> static void Exact(const std::vector<Value> &a, const std::vector<Value> &b)
{
	Require(a.size() == b.size() && (a.empty() ||
		std::memcmp(a.data(), b.data(), a.size() * sizeof(Value)) == 0),
		"initial state or metadata changed");
}

static void Exact(const CartesianSubgridFloat::TailInfo &a, const CartesianSubgridFloat::TailInfo &b)
{
	Require(a.terms == b.terms && a.l1 == b.l1 && a.relativeBound == b.relativeBound,
		"global initial-state FTZ ledger changed");
}

static void ExactMetadata(const CartesianSubgridDeviceData &a, const CartesianSubgridDeviceData &b)
{
	for (const auto &pair : {
		std::make_pair(&a.blockStart, &b.blockStart), {&a.blockEdge, &b.blockEdge},
		{&a.bandwidth, &b.bandwidth}, {&a.lowerStart, &b.lowerStart},
		{&a.sourceEdge, &b.sourceEdge}, {&a.portStart, &b.portStart}, {&a.portEdge, &b.portEdge},
		{&a.responseStart, &b.responseStart}, {&a.responseEdge, &b.responseEdge},
		{&a.sourceOffStart, &b.sourceOffStart}, {&a.sourceOffEdge, &b.sourceOffEdge},
		{&a.curlStart, &b.curlStart}, {&a.curlColumn, &b.curlColumn},
		{&a.transposeStart, &b.transposeStart}, {&a.transposeColumn, &b.transposeColumn},
		{&a.curlGeometry, &b.curlGeometry}})
	{
		Exact(*pair.first, *pair.second);
	}
	for (const auto &pair : {
		std::make_pair(&a.lower, &b.lower), {&a.sourceWeight, &b.sourceWeight},
		{&a.portWeight, &b.portWeight}, {&a.responseWeight, &b.responseWeight},
		{&a.rootConductance, &b.rootConductance}, {&a.portCurrent, &b.portCurrent},
		{&a.portLower, &b.portLower}, {&a.sourceOffWeight, &b.sourceOffWeight},
		{&a.ledger, &b.ledger}, {&a.curlWeight, &b.curlWeight}, {&a.curlLow, &b.curlLow},
		{&a.transposeWeight, &b.transposeWeight}, {&a.transposeLow, &b.transposeLow},
		{&a.curlMetrics, &b.curlMetrics}})
	{
		Exact(*pair.first, *pair.second);
	}
	Exact(a.stateTails, b.stateTails);
	Exact(a.factorTails, b.factorTails);
	Require(a.streamCurls == b.streamCurls && a.trackSourceWork == b.trackSourceWork &&
		a.stableTimestep == b.stableTimestep && a.portTimestep == b.portTimestep,
		"initial metadata changed source/port scheduling");
}

static Grid Make(unsigned int mode, bool ports, Grid::CurlLayout layout)
{
	const auto base = Base();
	const Grid::Index ratio{2, 4, 2};
	Grid grid(base, epsilon, mu, {{{{2, 2, 2}, {3, 4, 4}}, ratio}}, false,
		{}, {}, {}, {}, layout, mode >= 2 ? Children(base, ratio) : Grid::Axes{});
	grid.SetElectric([&](unsigned int n, const Grid::Point &p)
	{
		return mode % 2 && p[0] > base[0][4] ? 1e-100 :
			std::sin(731 * p[0] + 423 * p[1] - 511 * p[2] + n);
	});
	grid.SetMagnetic([&](unsigned int n, const Grid::Point &p)
	{
		return mode % 2 && p[1] > base[1][4] ? -1e-100 :
			std::cos(933 * p[0] - 331 * p[1] + 211 * p[2] + n) / 376.730313668;
	});
	if (ports)
	{
		grid.SetResistors({{{base[0][1], base[1][1], base[2][1]},
			{base[0][1], base[1][1], base[2][2]}, 50, .001}}, .8 * grid.StableTimestep());
	}
	grid.SetSources({{grid.EdgeCount() / 2, .1}});
	return grid;
}

static void Control(unsigned int mode, bool ports, Grid::CurlLayout layout)
{
	auto grid = Make(mode, ports, layout);
	const auto full = grid.DeviceStateData();
	const auto metadata = grid.DeviceSetupMetadata();
	Require(metadata.electric.empty() && metadata.magnetic.empty() &&
		metadata.electric.capacity() == 0 && metadata.magnetic.capacity() == 0,
		"setup metadata retained a global initial-state allocation");
	ExactMetadata(full, metadata);
	Require((metadata.stateTails.terms != 0) == (mode % 2 != 0),
		"initial metadata fixture did not exercise the intended global tail bound");
	const std::vector<double> cuts{Base()[0][3], Base()[0][5]};
	const Plan streamed(grid, 0, cuts), oracle(grid, 0, cuts, Plan::CurlPreparation::ExpandedOracle);
	Exact(streamed.StateTails(), full.stateTails);
	Exact(streamed.StateTails(), oracle.StateTails());
	for (std::size_t worker = 0; worker < streamed.Domains().size(); ++worker)
	{
		const auto &a = streamed.Domains()[worker], &b = oracle.Domains()[worker];
		Exact(a.globalElectric, b.globalElectric);
		Exact(a.globalMagnetic, b.globalMagnetic);
		Exact(a.fields.electric, b.fields.electric);
		Exact(a.fields.magnetic, b.fields.magnetic);
		for (std::size_t k = 0; k < a.globalElectric.size(); ++k)
		{
			Require(a.fields.electric[k] == full.electric[a.globalElectric[k]],
				"per-domain electric conversion changed the global initial state");
		}
		for (std::size_t k = 0; k < a.globalMagnetic.size(); ++k)
		{
			Require(a.fields.magnetic[k] == full.magnetic[a.globalMagnetic[k]],
				"per-domain magnetic conversion changed the global initial state");
		}
	}
	const auto after = grid.DeviceStateData();
	Exact(after.electric, full.electric);
	Exact(after.magnetic, full.magnetic);
	ExactMetadata(after, full);
	std::cout << "INITIAL_METADATA mode=" << mode << " ports=" << ports
		<< " spatial=" << (layout == Grid::CurlLayout::SpatialBlocks)
		<< " exact_metadata=1 empty_fields=1 exact_domain_state=1 source_unchanged=1" << std::endl;
}

static void Guards()
{
	using namespace CartesianSubgridFloat;
	TailInfo stored, scanned;
	const std::vector<double> values{1, -0., -.25, 1e-100, -1e-41, 0};
	const auto state = State(values, stored);
	std::vector<float> collected;
	StateTo(values, scanned, [&](float value) { collected.push_back(value); });
	Exact(state, std::vector<float>{1, -0.f, -.25f, 0, 0, 0});
	Exact(state, collected);
	Exact(stored, scanned);
	Require(scanned.terms == 2 && scanned.l1 == 1e-100 + 1e-41 && scanned.relativeBound == scanned.l1,
		"state scan changed the explicit negligible-tail oracle");
	for (const auto &bad : {std::vector<double>{1e-40}, {1e-25, 1e-39, 1e-39}, {1e50},
		{std::numeric_limits<double>::quiet_NaN()}, {std::numeric_limits<double>::infinity()}})
	{
		TailInfo info;
		Reject([&] { StateTo(bad, info, [](float) {}); });
		Require(info.terms == 0 && info.l1 == 0 && info.relativeBound == 0,
			"rejected state scan published a partial FTZ ledger");
	}
	auto grid = Make(0, false, Grid::CurlLayout::Patterns);
	grid.UpdateElectric(.1 * grid.StableTimestep());
	Reject([&] { grid.DeviceSetupMetadata(); });
	grid.UpdateMagnetic(.1 * grid.StableTimestep());
	grid.SetElectric([](unsigned int, const Grid::Point &) { return 1e-50; });
	Reject([&] { grid.DeviceStateData(); });
	Reject([&] { grid.DeviceSetupMetadata(); });
	Reject([&] { Plan invalid(grid, 0, {Base()[0][4]}); });
	Reject([&] { Plan invalid(grid, 0, {Base()[0][4]}, Plan::CurlPreparation::ExpandedOracle); });
	std::cout << "INITIAL_METADATA_GUARDS scalar_cases=5 partial_epoch=1 global_underflow=1" << std::endl;
}

int main()
{
	try
	{
		Guards();
		for (auto layout : {Grid::CurlLayout::Patterns, Grid::CurlLayout::SpatialBlocks})
		{
			for (unsigned int mode = 0; mode < 4; ++mode)
			{
				for (bool ports : {false, true})
				{
					Control(mode, ports, layout);
				}
			}
		}
		std::cout << "PASS: metadata-only global state checks and exact per-domain initial fields" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

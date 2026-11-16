/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
using Plan = CartesianSubgridDomainPlan;

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

template <class T> static void Exact(const std::vector<T> &a, const std::vector<T> &b)
{
	Require(a.size() == b.size() &&
		(a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0),
		"streamed domain preparation differs from the expanded oracle");
}

static void Exact(const CartesianSubgridDeviceData &a, const CartesianSubgridDeviceData &b)
{
	for (const auto &pair : {
		std::make_pair(&a.curlStart, &b.curlStart), {&a.curlColumn, &b.curlColumn},
		{&a.transposeStart, &b.transposeStart}, {&a.transposeColumn, &b.transposeColumn},
		{&a.blockStart, &b.blockStart}, {&a.blockEdge, &b.blockEdge}, {&a.bandwidth, &b.bandwidth},
		{&a.lowerStart, &b.lowerStart}, {&a.sourceEdge, &b.sourceEdge}, {&a.portStart, &b.portStart},
		{&a.portEdge, &b.portEdge}, {&a.responseStart, &b.responseStart}, {&a.responseEdge, &b.responseEdge},
		{&a.sourceOffStart, &b.sourceOffStart}, {&a.sourceOffEdge, &b.sourceOffEdge},
		{&a.curlGeometry, &b.curlGeometry}})
	{
		Exact(*pair.first, *pair.second);
	}
	for (const auto &pair : {
		std::make_pair(&a.electric, &b.electric), {&a.magnetic, &b.magnetic},
		{&a.curlWeight, &b.curlWeight}, {&a.curlLow, &b.curlLow},
		{&a.transposeWeight, &b.transposeWeight}, {&a.transposeLow, &b.transposeLow},
		{&a.lower, &b.lower}, {&a.sourceWeight, &b.sourceWeight}, {&a.portWeight, &b.portWeight},
		{&a.responseWeight, &b.responseWeight}, {&a.rootConductance, &b.rootConductance},
		{&a.portCurrent, &b.portCurrent}, {&a.portLower, &b.portLower},
		{&a.sourceOffWeight, &b.sourceOffWeight}, {&a.ledger, &b.ledger}, {&a.curlMetrics, &b.curlMetrics}})
	{
		Exact(*pair.first, *pair.second);
	}
	Require(a.streamCurls == b.streamCurls && a.trackSourceWork == b.trackSourceWork &&
		a.stableTimestep == b.stableTimestep && a.portTimestep == b.portTimestep,
		"streamed domain changed field/port scheduling");
}

static Grid Fixture(bool multiple)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		axes[n].push_back(0);
		for (unsigned int i = 0; i < 12; ++i)
		{
			axes[n].push_back(axes[n].back() + .001 * (1 + .1 * n + .01 * i));
		}
	}
	const std::vector<Grid::Region> regions = multiple ?
		std::vector<Grid::Region>{{{{2, 2, 2}, {5, 9, 9}}, {2, 2, 1}},
			{{{7, 3, 3}, {10, 10, 10}}, {1, 2, 4}}} :
		std::vector<Grid::Region>{{{{3, 3, 3}, {9, 9, 9}}, {2, 4, 2}}};
	Grid grid(axes, 8.85418781762e-12, 1.2566370614359173e-6, regions, false,
		[&](const Grid::Point &, const Grid::Point &stop)
		{
			return Grid::Medium{8.85418781762e-12 * (stop[0] <= axes[0][6] ? 2.25 : 1),
				1.2566370614359173e-6};
		},
		[&](unsigned int, const Grid::Point &start, const Grid::Point &stop)
		{
			return start[2] == axes[2][2] && stop[2] == axes[2][2];
		}, Grid::BoundaryConditions{}, {}, Grid::CurlLayout::SpatialBlocks);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return std::cos(933 * p[0] - 1731 * p[1] + 411 * p[2] + n) / 376.730313668;
	});
	grid.SetResistors({{{axes[0][1], axes[1][4], axes[2][4]},
		{axes[0][1], axes[1][4], axes[2][5]}, 50, .001},
		{{axes[0][10], axes[1][6], axes[2][6]},
		{axes[0][10], axes[1][6], axes[2][7]}, 75, 0}}, .8 * grid.StableTimestep());
	grid.SetSources({{grid.EdgeCount() / 2, .1}, {grid.EdgeCount() / 3, -.03}});
	return grid;
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool multiple)
{
	auto grid = Fixture(multiple), first = grid, second = grid;
	const double span = grid.BaseLines()[axis].back();
	const std::vector<double> cuts{.31 * span, .68 * span};
	auto streamed = std::make_shared<Plan>(grid, axis, cuts);
	auto expanded = std::make_shared<Plan>(grid, axis, cuts, Plan::CurlPreparation::ExpandedOracle);
	Require(streamed->ElectricCount() == expanded->ElectricCount() &&
		streamed->MagneticCount() == expanded->MagneticCount() &&
		streamed->ReassignedElectricCount() == expanded->ReassignedElectricCount(),
		"streamed preparation changed ownership counts");
	for (std::size_t k = 0; k < streamed->Domains().size(); ++k)
	{
		const auto &a = streamed->Domains()[k], &b = expanded->Domains()[k];
		Require(a.ownedElectric == b.ownedElectric && a.ownedMagnetic == b.ownedMagnetic,
			"streamed preparation changed owned fields");
		Exact(a.globalElectric, b.globalElectric);
		Exact(a.globalMagnetic, b.globalMagnetic);
		Exact(a.sendElectric, b.sendElectric);
		Exact(a.sendMagnetic, b.sendMagnetic);
		Exact(a.fields, b.fields);
	}
	for (auto field : {CartesianSubgridField::Electric, CartesianSubgridField::Magnetic})
	{
		const auto &a = streamed->Transfers(field), &b = expanded->Transfers(field);
		Require(a.size() == b.size(), "streamed preparation changed exchange counts");
		for (std::size_t k = 0; k < a.size(); ++k)
		{
			Require(a[k].sender == b[k].sender && a[k].receiver == b[k].receiver &&
				a[k].sendOffset == b[k].sendOffset && a[k].receiveOffset == b[k].receiveOffset,
				"streamed preparation changed a halo transfer");
		}
	}
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend *> pointers;
	for (std::size_t k = 0; k < streamed->Domains().size(); ++k)
	{
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "requested domain-stream backend unavailable");
		target.Verify(*backends.back());
		pointers.push_back(backends.back().get());
	}
	auto candidate = CreateCartesianSubgridDomains(streamed, pointers);
	auto oracle = CreateCartesianSubgridDomains(expanded, pointers);
	const double dt = .8 * grid.StableTimestep();
	for (unsigned int k = 0; k < 256; ++k)
	{
		const double amplitude = std::sin(.17 * k);
		candidate->Step(dt, amplitude);
		oracle->Step(dt, amplitude);
	}
	candidate->ReadFields(first);
	oracle->ReadFields(second);
	const auto a = first.DeviceStateData(), b = second.DeviceStateData();
	Exact(a.electric, b.electric);
	Exact(a.magnetic, b.magnetic);
	std::cout << "DOMAIN_STREAM axis=" << axis << " multiple=" << multiple
		<< " leaves=" << grid.CellCount() << " metadata_exact=1 fields_exact=1" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SubgridDeviceTarget target;
		if (argc == 2)
		{
			target = SubgridDeviceTarget::FromDomainFlag(argv[1]);
		}
		Require(argc <= 2, "expected at most one backend flag");
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			for (bool multiple : {false, true})
			{
				Control(target, axis, multiple);
			}
		}
		std::cout << "PASS: streamed domain metadata and fields match the expanded oracle" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

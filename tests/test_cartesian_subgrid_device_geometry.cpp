/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/subgrid_device_geometry.h"
#include "subgrid_device_target.h"
#include "subgrid_fp32_oracle.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>

using Grid = CartesianSubgrid;
static constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}
static double Error(const std::vector<float> &a, const std::vector<float> &b)
{
	Require(a.size() == b.size(), "field dimensions differ");
	double delta = 0, scale = 0;
	for (std::size_t i = 0; i < a.size(); ++i)
	{
		Require(std::isfinite(b[i]), "nonfinite structured device field");
		delta = std::max(delta, std::abs(double(a[i]) - b[i]));
		scale = std::max(scale, std::abs(double(a[i])));
	}
	return delta / std::max(scale, 1e-300);
}
static Grid Fixture(const std::string &name, Grid::CurlLayout layout)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i <= 10; ++i)
		{
			axes[n].push_back(.001 * i * (1 + .1 * n) + (name == "graded" ? 1e-6 * i * i : 0));
		}
	}
	std::vector<Grid::Region> regions{{{{3, 3, 3}, {7, 7, 7}}, name == "uniform" ? Grid::Index{1, 1, 1} : Grid::Index{2, 2, 2}}};
	if (name == "multiple")
	{
		regions = {{{{1, 1, 1}, {3, 3, 3}}, {2, 4, 1}}, {{{6, 6, 6}, {9, 9, 9}}, {4, 1, 2}}};
	}
	Grid::MediumSampler material;
	if (name == "material")
	{
		material = [&](const Grid::Point &lo, const Grid::Point &)
		{
			return Grid::Medium{epsilon * (lo[0] >= axes[0][5] ? 2.3 : 1), mu};
		};
	}
	const auto inside = [&](const Grid::Point &lo, const Grid::Point &hi)
	{
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (lo[n] < axes[n][4] - 1e-14 || hi[n] > axes[n][6] + 1e-14)
			{
				return false;
			}
		}
		return true;
	};
	Grid::OpaqueSelector opaque;
	if (name == "opaque")
	{
		opaque = inside;
	}
	Grid::PECSelector pec;
	if (name == "pec")
	{
		pec = [&](unsigned int, const Grid::Point &a, const Grid::Point &b) { return inside(a, b); };
	}
	Grid grid(axes, epsilon, mu, regions, false, material, pec, {}, opaque, layout);
	grid.SetElectric([](unsigned int n, const Grid::Point &p) { return std::sin(731 * p[(n + 1) % 3]) + .2 * std::cos(451 * p[(n + 2) % 3]); });
	if (name != "opaque")
	{
		grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return .001 * std::cos(317 * p[n]); });
	}
	return grid;
}

template <class Action> static void Reject(Action action)
{
	try
	{
		action();
	}
	catch (const std::invalid_argument &)
	{
		return;
	}
	throw std::runtime_error("malformed spatial metadata or arithmetic was admitted");
}
static void NegativeControls()
{
	const auto grid = Fixture("graded", Grid::CurlLayout::SpatialBlocks);
	const auto data = grid.DeviceData();
	const auto check = [&](const CartesianSubgridDeviceData &invalid)
	{
		CartesianSubgridDeviceGeometry::Validate(invalid.curlGeometry, invalid.curlMetrics,
			invalid.electric.size(), invalid.magnetic.size());
		CartesianSubgridDeviceGeometry::ValidateRows(invalid);
	};
	check(data);
	auto invalid = data;
	invalid.curlGeometry[0] = 0;
	Reject([&] { check(invalid); });
	invalid = data;
	invalid.curlGeometry[3] += 1;
	Reject([&] { check(invalid); });
	invalid = data;
	invalid.curlGeometry[4 + 9 * invalid.curlGeometry[0] + 3] += 1;
	Reject([&] { check(invalid); });
	invalid = data;
	invalid.curlMetrics[0] = std::numeric_limits<float>::quiet_NaN();
	Reject([&] { check(invalid); });
	invalid = data;
	invalid.curlMetrics[0] = std::nextafter(invalid.curlMetrics[0], std::numeric_limits<float>::infinity());
	Reject([&] { check(invalid); });
	invalid = data;
	bool changed = false;
	for (std::size_t face = 0; face < invalid.magnetic.size() && !changed; ++face)
	{
		if (subgrid_geometry_regular(invalid.curlGeometry.data(), std::uint32_t(face)))
		{
			const auto edge = subgrid_geometry_edge_index(invalid.curlGeometry.data(), std::uint32_t(face), 0);
			for (std::size_t k = invalid.transposeStart[edge]; k < invalid.transposeStart[edge + 1]; ++k)
			{
				if (invalid.transposeColumn[k] == face)
				{
					invalid.transposeWeight[k] = std::nextafter(invalid.transposeWeight[k], std::numeric_limits<float>::infinity());
					changed = true;
					break;
				}
			}
		}
	}
	Require(changed, "negative transpose control found no structured face");
	Reject([&] { check(invalid); });
	invalid = data;
	changed = false;
	for (std::size_t face = 0; face < invalid.magnetic.size() && !changed; ++face)
	{
		if (subgrid_geometry_regular(invalid.curlGeometry.data(), std::uint32_t(face)))
		{
			auto &descriptor = invalid.curlGeometry[invalid.curlGeometry[3] + face];
			descriptor = (descriptor & ~12u) | ((descriptor & 3u) << 2);
			changed = true;
		}
	}
	Require(changed, "negative ordering control found no structured face");
	Reject([&] { check(invalid); });
	std::cout << "PASS: malformed spatial headers, prefixes, coefficients and transpose pairs rejected" << std::endl;
}
// Independent expanded r06 preparation; intentionally retains CSR in this test.
static CartesianSubgridDeviceOperators ExpandedOperators(const CartesianSubgridDeviceData &data)
{
	CartesianSubgridDeviceOperators result;
	CartesianSubgridDeviceCurl transpose(data.transposeStart, data.transposeColumn,
		data.transposeWeight, data.transposeLow, data.magnetic.size());
	CartesianSubgridDeviceCurl curl(data.curlStart, data.curlColumn,
		data.curlWeight, data.curlLow, data.electric.size());
	const auto spatialRows = CartesianSubgridDeviceGeometry::Validate(
		data.curlGeometry, data.curlMetrics, data.electric.size(), data.magnetic.size());
	if (spatialRows)
	{
		CartesianSubgridDeviceGeometry::ValidateRows(data);
		std::vector<std::uint32_t> start{0}, columns;
		std::vector<float> high, low;
		for (std::size_t row = 0; row < data.magnetic.size(); ++row)
		{
			if (!subgrid_geometry_regular(data.curlGeometry.data(), std::uint32_t(row)))
			{
				for (std::size_t k = data.curlStart[row]; k < data.curlStart[row + 1]; ++k)
				{
					columns.push_back(data.curlColumn[k]);
					high.push_back(data.curlWeight[k]);
					low.push_back(data.curlLow[k]);
				}
				start.push_back(static_cast<std::uint32_t>(columns.size()));
			}
		}
		CartesianSubgridDeviceCurl packedCurl(start, columns, high, low, data.electric.size());
		high = data.transposeWeight;
		low = data.transposeLow;
		for (std::size_t edge = 0; edge < data.electric.size(); ++edge)
		{
			for (std::size_t k = data.transposeStart[edge]; k < data.transposeStart[edge + 1]; ++k)
			{
				const auto face = data.transposeColumn[k];
				if (!subgrid_geometry_regular(data.curlGeometry.data(), face))
				{
					continue;
				}
				// A regular transpose term stores an exact 0..3 slot, not an
				// approximate coefficient. Both gathers read the same pair.
				unsigned int term = 0;
				while (term < 4 && subgrid_geometry_edge_index(data.curlGeometry.data(), face, term) != edge)
				{
					++term;
				}
				if (term == 4)
				{
					throw std::logic_error("CartesianSubgrid: validated transpose lost its spatial edge slot");
				}
				high[k] = static_cast<float>(term);
				low[k] = 0;
			}
		}
		CartesianSubgridDeviceCurl packedTranspose(data.transposeStart, data.transposeColumn,
			high, low, data.magnetic.size());
		const std::size_t oldBytes = curl.DeviceBytes() + transpose.DeviceBytes();
		const std::size_t packedBytes = packedCurl.DeviceBytes() + packedTranspose.DeviceBytes() +
			data.curlGeometry.size() * sizeof(std::uint32_t) + data.curlMetrics.size() * sizeof(float);
		// A simple uniform grid may already have smaller translated patterns.
		if (packedBytes < oldBytes)
		{
			curl = std::move(packedCurl);
			transpose = std::move(packedTranspose);
			result.spatialRows = spatialRows;
			result.legacyBytes = oldBytes;
			result.geometry = data.curlGeometry;
			result.metrics = data.curlMetrics;
		}
	}

	result.curl = std::move(curl);
	result.transpose = std::move(transpose);
	return result;
}

template <class T> static bool Exact(const std::vector<T> &a, const std::vector<T> &b)
{
	return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0);
}
static void CompareCurl(const CartesianSubgridDeviceCurl &a, const CartesianSubgridDeviceCurl &b)
{
	Require(a.start == b.start && a.offset == b.offset && a.base == b.base && a.pattern == b.pattern &&
		Exact(a.weight, b.weight) && Exact(a.residual, b.residual),
		"streamed operator differs from the expanded preparation");
}
static void CheckStream(const Grid &grid)
{
	const auto data = grid.DeviceData();
	const auto expected = ExpandedOperators(data);
	const auto actual = grid.DeviceCurls();
	CompareCurl(expected.curl, actual.curl);
	CompareCurl(expected.transpose, actual.transpose);
	Require(expected.geometry == actual.geometry && Exact(expected.metrics, actual.metrics) &&
		expected.spatialRows == actual.spatialRows,
		"streamed spatial pool or byte admission changed");
	const auto state = grid.DeviceStateData();
	Require(state.streamCurls && !data.streamCurls && state.curlColumn.empty() &&
		state.curlStart == std::vector<std::uint32_t>{0} && state.transposeStart.empty() &&
		state.transposeColumn.empty() && state.curlWeight.empty() && state.curlLow.empty() &&
		state.transposeWeight.empty() && state.transposeLow.empty() &&
		state.curlGeometry.empty() && state.curlMetrics.empty(),
		"state-only preparation retained expanded curls");
	Require(Exact(data.electric, state.electric) && Exact(data.magnetic, state.magnetic) &&
		Exact(data.lower, state.lower) && Exact(data.sourceWeight, state.sourceWeight) &&
		Exact(data.portWeight, state.portWeight) && Exact(data.responseWeight, state.responseWeight) &&
		Exact(data.rootConductance, state.rootConductance) && Exact(data.portCurrent, state.portCurrent) &&
		Exact(data.portLower, state.portLower) && Exact(data.sourceOffWeight, state.sourceOffWeight) &&
		Exact(data.ledger, state.ledger) && data.blockStart == state.blockStart &&
		data.blockEdge == state.blockEdge && data.bandwidth == state.bandwidth &&
		data.lowerStart == state.lowerStart && data.sourceEdge == state.sourceEdge &&
		data.portStart == state.portStart && data.portEdge == state.portEdge &&
		data.responseStart == state.responseStart && data.responseEdge == state.responseEdge &&
		data.sourceOffStart == state.sourceOffStart && data.sourceOffEdge == state.sourceOffEdge &&
		data.stableTimestep == state.stableTimestep && data.portTimestep == state.portTimestep,
		"state-only preparation changed fields, factors, sources or ports");
}
static bool Control(const SubgridDeviceTarget &target, const std::string &name, unsigned int steps, bool source)
{
	auto cpu = Fixture(name, Grid::CurlLayout::Patterns);
	auto oracle = cpu;
	auto device = Fixture(name, Grid::CurlLayout::SpatialBlocks);
	if (source)
	{
		for (Grid *grid : {&cpu, &oracle, &device})
		{
			grid->SetSources({{0, .1}});
		}
	}
	CheckStream(cpu);
	CheckStream(device);
	const auto expected = cpu.DeviceData(), prepared = device.DeviceData();
	Require(expected.curlWeight == prepared.curlWeight && expected.curlLow == prepared.curlLow &&
		expected.curlColumn == prepared.curlColumn, "device geometry changed the expanded oracle");
	Require(CartesianSubgridDeviceGeometry::Validate(prepared.curlGeometry, prepared.curlMetrics,
		prepared.electric.size(), prepared.magnetic.size()) > 0, "spatial packing has no admitted rows");
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "selected physical backend unavailable; no fallback");
	target.Verify(*backend);
	auto execution = backend->CreateCartesianSubgrid(device);
	auto reference = CreateCartesianSubgridReference(oracle);
	const double dt = .8 * cpu.StableTimestep(), energy = cpu.Energy(dt);
	for (unsigned int step = 0; step < steps; ++step)
	{
		const double drive = source ? .2 * std::exp(-std::pow((double(step) - 20) / 8, 2)) : 0;
		cpu.Step(dt, drive);
		execution->Step(dt, drive);
		reference->Step(dt, drive);
	}
	execution->ReadFields(device);
	reference->ReadFields(oracle);
	const auto a = cpu.DeviceData(), b = device.DeviceData(), c = oracle.DeviceData();
	const double error = std::max(Error(a.electric, b.electric), Error(a.magnetic, b.magnetic));
	const double oracleError = std::max(Error(c.electric, b.electric), Error(c.magnetic, b.magnetic));
	const double drift = std::abs(device.Energy(dt) - energy) / energy;
	std::cout << "DEVICE_GEOMETRY " << name << " steps=" << steps << " CPU_error=" << error
		<< " FP32_error=" << oracleError << " drift=" << drift << " structured_rows=" << execution->StructuredRows()
		<< " device_bytes=" << execution->StorageBytes() << std::endl;
	Require(error < 2e-4, "structured device/CPU field error exceeded the unchanged 2e-4 gate");
	Require(oracleError < 2e-6, "structured device/FP32 oracle exceeded the unchanged 2e-6 gate");
	if (!source)
	{
		Require(drift < 5e-4, "structured device invariant drift exceeded 5e-4");
	}
	return execution->StructuredRows() != 0;
}
int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "usage: test_cartesian_subgrid_device_geometry [--metal|--cuda]");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromFlag(argv[1]) : SubgridDeviceTarget{};
		NegativeControls();
		std::size_t used = 0;
		for (const std::string name : {"uniform", "local", "graded", "multiple", "material", "pec", "opaque"})
		{
			used += Control(target, name, 2000, true);
		}
		used += Control(target, "graded", 20000, false);
		Require(used >= 2, "all structured kernels were bypassed by the memory admission gate");
		std::cout << "PASS: shared geometric curl/transpose, exact sparse exceptions and actual structured device fields" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

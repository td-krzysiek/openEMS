/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "FDTD/gpu_cartesiansubgrid_dft.h"
#include "FDTD/kernels/subgrid_dft.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

using Grid = CartesianSubgrid;
using Complex = std::complex<double>;
static void Require(bool condition, const char *message)
{
	if (!condition)
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
	throw std::runtime_error("invalid DFT operation accepted");
}

//! Cancellation oracle does not use the DFT helper or device execution.
static void LongSum()
{
	float field = .1f, sums[4] = {};
	const std::uint32_t index = 0;
	Op_SubgridDFT op{&field, &index, sums, 1, 0, 0, 0};
	for (unsigned int k = 0; k < 99999; ++k)
	{
		op.real = k % 3 == 0 ? 1e8f : k % 3 == 1 ? 1 : -1e8f;
		op.imag = -op.real;
		op(0);
	}
	const double expected = double(field) * 33333;
	Require(std::abs(double(sums[0]) + sums[1] - expected) < 1e-6 * expected &&
		std::abs(double(sums[2]) + sums[3] + expected) < 1e-6 * expected,
		"long compensated DFT cancellation lost small terms");
	std::cout << "DFT_LONG samples=99999 cancellation=1\n";
}

static void Helper(const SubgridDeviceTarget &target)
{
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "DFT helper backend missing");
	target.Verify(*backend);
	backend->Activate();
	const std::vector<float> input{.1f, -.3f, 1e-18f, 1e18f};
	auto *field = static_cast<float *>(backend->DeviceAlloc(input.size() * sizeof(float), input.data()));
	const std::vector<std::uint32_t> indices{3, 0, 1, 0, 2};
	Reject([&] { CartesianSubgridDFTBytes(1, 0); });
	Reject([&] { CartesianSubgridDFTBytes(std::numeric_limits<std::uint32_t>::max(), std::numeric_limits<std::uint32_t>::max()); });
	Reject([&] { GPU_CreateCartesianSubgridDFT(*backend, field, input.size(), {4}, 2); });
	auto dft = GPU_CreateCartesianSubgridDFT(*backend, field, input.size(), indices, 2);
	auto empty = GPU_CreateCartesianSubgridDFT(*backend, nullptr, 0, {}, 2);
	Require(dft->StorageBytes() == indices.size() * 36 && empty->StorageBytes() == 0,
		"DFT payload/empty allocation accounting differs");
	std::vector<Complex> sums;
	dft->Read(0, sums);
	Require(std::all_of(sums.begin(), sums.end(), [](Complex value) { return value == Complex{}; }), "nonzero initial DFT");
	Reject([&] { dft->Accumulate({{1, 0}}); });
	Reject([&] { dft->Accumulate({{std::numeric_limits<float>::infinity(), 0}, {0, 0}}); });
	Reject([&] { dft->Read(2, sums); });
	const auto before = dft->DownloadBytes();
	for (unsigned int k = 0; k < 1023; ++k)
	{
		const float weight = k % 3 == 0 ? 1e8f : k % 3 == 1 ? 1 : -1e8f;
		dft->Accumulate({{weight, -weight}, {-weight, weight}});
		empty->Accumulate({{weight, -weight}, {-weight, weight}});
	}
	Require(dft->DownloadBytes() == before && empty->DownloadBytes() == 0,
		"DFT accumulation downloaded data");
	for (std::size_t frequency = 0; frequency < 2; ++frequency)
	{
		dft->Read(frequency, sums);
		for (std::size_t k = 0; k < indices.size(); ++k)
		{
			const double value = double(input[indices[k]]) * 341 * (frequency ? -1 : 1);
			Require(std::abs(sums[k] - Complex(value, -value)) < 2e-6 * std::abs(value),
				"device DFT cancellation, duplicate order or scaling differs");
		}
	}
	empty->Read(1, sums);
	Require(sums.empty() && empty->DownloadBytes() == 0, "empty DFT read allocated/transferred values");
	Require(dft->DownloadBytes() - before == 2 * indices.size() * 16, "DFT sum download accounting differs");
	std::cout << "DFT_HELPER samples=1023 scales=4 duplicate=1 no_sample_download=1 backend=" << backend->GetName() << '\n';
}

static Grid Fixture(bool ports)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int k = 0; k <= 8; ++k)
		{
			axes[n].push_back(.001 * k * (1 + .03 * n + .01 * k));
		}
	}
	Grid grid(axes, 8.85418781762e-12, 1.2566370614359173e-6,
		std::vector<Grid::Region>{{{{2, 2, 2}, {6, 6, 6}}, {2, 4, 2}}}, false);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(791 * p[0] - 1183 * p[1] + 913 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return .001 * std::cos(517 * p[0] + 1233 * p[1] - 723 * p[2] + n);
	});
	if (ports)
	{
		grid.SetResistors({{{axes[0][1], axes[1][1], axes[2][1]},
			{axes[0][1], axes[1][1], axes[2][2]}, 50, 0}}, .7 * grid.StableTimestep());
	}
	grid.SetSources({{grid.EdgeCount() / 2, .1}});
	return grid;
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool ports)
{
	auto original = Fixture(ports), one = original, split = original, oracleGrid = original;
	const double timestep = .7 * original.StableTimestep();
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "DFT single backend missing");
	target.Verify(*backend);
	auto single = backend->CreateCartesianSubgrid(one);
	auto oracle = backend->CreateCartesianSubgrid(oracleGrid);
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<GPU_Backend *> pointers;
	for (unsigned int k = 0; k < 3; ++k)
	{
		target.SelectDomain(k);
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "DFT domain backend missing");
		target.Verify(*backends.back());
		pointers.push_back(backends.back().get());
	}
	const auto &lines = original.BaseLines()[axis];
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(split, axis, std::vector<double>{lines[3], lines[5]});
	auto distributed = CreateCartesianSubgridDomains(plan, pointers);
	const std::vector<std::uint32_t> selection[2] = {
		{static_cast<std::uint32_t>(original.EdgeCount() - 1), 0, 13, 13, static_cast<std::uint32_t>(original.EdgeCount() / 2)},
		{static_cast<std::uint32_t>(original.FaceCount() - 1), 0, 11, 11, static_cast<std::uint32_t>(original.FaceCount() / 2)}};
	constexpr std::size_t budget = 1 << 20;
	std::vector<Complex> values;
	for (auto *execution : {single.get(), distributed.get()})
	{
		Reject([&] { execution->AccumulateObservationDFT(0, {{1, 0}}); });
		Reject([&] { execution->AddObservationDFT(true, selection[0], 3, 1); });
		Reject([&] { execution->AddObservationDFT(true, {static_cast<std::uint32_t>(original.EdgeCount())}, 3, budget); });
		Require(execution->AddObservationDFT(true, selection[0], 3, budget) == 0 &&
			execution->AddObservationDFT(false, selection[1], 3, budget) == 1 &&
			execution->AddObservationDFT(true, {}, 1, budget) == 2, "DFT group registration order differs");
		Reject([&] { execution->ReadObservationDFT(0, 3, values); });
		execution->ReadObservationDFT(0, 0, values);
		Require(values == std::vector<Complex>(selection[0].size()), "nonzero execution initial DFT");
		execution->ReadObservationDFT(2, 0, values);
		Require(values.empty(), "empty execution DFT nonempty");
	}
	const auto singleHost = single->HostObservationStorageBytes(), splitHost = distributed->HostObservationStorageBytes();
	std::vector<std::vector<Complex>> expected[2];
	for (unsigned int field = 0; field < 2; ++field)
	{
		expected[field].assign(3, std::vector<Complex>(selection[field].size()));
	}
	double error = 0;
	for (unsigned int step = 0; step < 128; ++step)
	{
		const double amplitude = .03 * std::sin(.17 * step);
		for (auto *execution : {single.get(), distributed.get()})
		{
			execution->UpdateElectric(timestep, amplitude);
			Reject([&] { execution->AccumulateObservationDFT(0, {{1, 0}, {0, 1}, {1, 1}}); });
			Reject([&] { execution->ReadObservationDFT(0, 0, values); });
			execution->UpdateMagnetic(timestep);
		}
		oracle->Step(timestep, amplitude);
		if (step % 3 != 0)
		{
			continue;
		}
		oracle->ReadFields(oracleGrid);
		const auto raw = oracleGrid.DeviceStateData();
		for (unsigned int field = 0; field < 2; ++field)
		{
			std::vector<std::complex<float>> weights;
			for (unsigned int frequency = 0; frequency < 3; ++frequency)
			{
				const double phase = (.013 + .007 * frequency) * (step + (field ? .5 : 0));
				weights.emplace_back(static_cast<float>(std::cos(phase)), static_cast<float>(-std::sin(phase)));
				for (std::size_t k = 0; k < selection[field].size(); ++k)
				{
					const auto &state = field ? raw.magnetic : raw.electric;
					expected[field][frequency][k] += Complex(weights.back()) * double(state[selection[field][k]]);
				}
			}
			for (auto *execution : {single.get(), distributed.get()})
			{
				const auto downloads = execution->ObservationDFTDownloadBytes();
				execution->AccumulateObservationDFT(field, weights);
				Require(execution->ObservationDFTDownloadBytes() == downloads, "DFT sample downloaded sums");
				Reject([&] { execution->AccumulateObservationDFT(field, weights); });
				Reject([&] { execution->AddObservationDFT(true, {}, 1, budget); });
			}
		}
		Require(single->HostObservationStorageBytes() == singleHost && distributed->HostObservationStorageBytes() == splitHost,
			"DFT sample materialized complete host fields");
		if (step % 24 == 0 || step == 126)
		{
			for (unsigned int field = 0; field < 2; ++field)
			{
				for (unsigned int frequency = 0; frequency < 3; ++frequency)
				{
					double scale = 0;
					for (const auto value : expected[field][frequency])
					{
						scale = std::max(scale, std::abs(value));
					}
					for (auto *execution : {single.get(), distributed.get()})
					{
						execution->ReadObservationDFT(field, frequency, values);
						Require(values.size() == selection[field].size(), "DFT selection count changed");
						for (std::size_t k = 0; k < values.size(); ++k)
						{
							error = std::max(error, std::abs(values[k] - expected[field][frequency][k]) / scale);
						}
					}
				}
			}
		}
	}
	Require(error < 2e-5, "DFT complex sum differs from independent eager-field accumulation");
	oracle->ReadFields(oracleGrid);
	single->ReadFields(one);
	distributed->ReadFields(split);
	const auto a = one.DeviceStateData(), b = split.DeviceStateData(), c = oracleGrid.DeviceStateData();
	Require(a.electric == c.electric && a.magnetic == c.magnetic, "DFT observations changed single-device physics");
	double fieldError = 0, fieldScale = 0;
	for (std::size_t k = 0; k < a.electric.size(); ++k)
	{
		fieldError = std::max(fieldError, double(std::abs(a.electric[k] - b.electric[k])));
		fieldScale = std::max(fieldScale, double(std::abs(a.electric[k])));
	}
	Require(fieldError < 2e-5 * fieldScale, "DFT domain physics differs from single control");
	std::cout << "DFT_CONTROL axis=" << axis << " ports=" << ports << " complex_error=" << error
		<< " no_sample_download=1 exact_single_physics=1\n";
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "unexpected DFT arguments");
		const auto target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		LongSum();
		Helper(target);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, false);
			Control(target, axis, true);
		}
		std::cout << "PASS: owned-state device DFT, bounded storage, complete epochs and domains\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

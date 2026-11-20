/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Grid = CartesianSubgrid;
static void Require(bool condition, const char *message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}
static void Control(double thickness, const SubgridDeviceTarget *target)
{
	Grid::Axes axes;
	for (unsigned int n = 0; n < 2; ++n)
	{
		for (unsigned int i = 0; i <= 4; ++i)
		{
			axes[n].push_back(.001 * i);
		}
	}
	axes[2] = {0, .0005, .0008, .0008 + thickness, .0015, .002};
	const std::vector<Grid::Region> regions{{{{1, 1, 1}, {3, 3, 4}}, {2, 2, 1}}};
	const auto opaque = [&](const Grid::Point &lo, const Grid::Point &hi)
	{
		return lo[0] >= .0015 && hi[0] <= .0025 && lo[1] >= .0015 && hi[1] <= .0025 &&
			lo[2] >= axes[2][2] && hi[2] <= axes[2][3];
	};
	constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;
	bool rejected = false;
	try
	{
		Grid cut(axes, epsilon, mu, {{{{1, 1, 1}, {3, 3, 3}}, {2, 2, 1}}}, false, {}, {}, {}, opaque);
	}
	catch (const std::invalid_argument &error)
	{
		rejected = std::string(error.what()).find("coarse interface trace") != std::string::npos;
	}
	Require(rejected, "unsupported conductor/coarse-fine boundary pattern was admitted");
	// An installed classifier with no excluded leaves must leave the old scheme exact.
	Grid unchanged(axes, epsilon, mu, regions), emptyMask(axes, epsilon, mu, regions,
		false, {}, {}, {}, [](const Grid::Point &, const Grid::Point &) { return false; });
	const auto seed = [](unsigned int n, const Grid::Point &p)
	{
		return std::sin(731 * p[0] - 433 * p[1] + 211 * p[2] + n);
	};
	unchanged.SetElectric(seed);
	emptyMask.SetElectric(seed);
	const double emptyStep = .4 * unchanged.StableTimestep();
	Require(emptyMask.StableTimestep() == unchanged.StableTimestep(), "empty volume mask changed timestep");
	for (unsigned int step = 0; step < 64; ++step)
	{
		unchanged.Step(emptyStep);
		emptyMask.Step(emptyStep);
	}
	const auto expectedData = unchanged.DeviceData(), emptyData = emptyMask.DeviceData();
	Require(expectedData.electric == emptyData.electric && expectedData.magnetic == emptyData.magnetic &&
		unchanged.Energy(emptyStep) == emptyMask.Energy(emptyStep) && emptyMask.OpaqueBoundaryFaces().empty(),
		"empty volume classifier changed unchanged fields or energy");
	Grid grid(axes, epsilon, mu, regions, false, {}, {}, {}, opaque);
	std::size_t opaqueCount = 0, surfaceEdges = 0;
	double excludedVolume = 0;
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		if (grid.LeafIsOpaque(leaf))
		{
			++opaqueCount;
			const auto length = grid.LeafLengths(leaf);
			excludedVolume += length[0] * length[1] * length[2];
			for (unsigned int n = 0; n < 3; ++n)
			{
				for (unsigned int i = 0; i < 2; ++i)
				{
					for (unsigned int j = 0; j < 2; ++j)
					{
						surfaceEdges += !grid.LeafEdgeVoltage(leaf, n, i, j).weights.empty();
					}
				}
			}
		}
	}
	const auto data = grid.DeviceData();
	auto invalidMagnetic = data.magnetic;
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		if (grid.LeafIsOpaque(leaf))
		{
			invalidMagnetic[6 * leaf] = 1;
			break;
		}
	}
	rejected = false;
	try
	{
		grid.ImportDeviceFields(data.electric.data(), data.electric.size(), invalidMagnetic.data(), invalidMagnetic.size());
	}
	catch (const std::invalid_argument &error)
	{
		rejected = std::string(error.what()).find("excluded conductor volume") != std::string::npos;
	}
	Require(rejected, "device import admitted an interior conductor field");
	Require(opaqueCount == 4 && surfaceEdges > 0 && grid.PECTraceCount() > 0,
		"opaque volume did not retain boundary traces and eliminate its interior");
	Require(std::abs(excludedVolume / (1e-6 * thickness) - 1) < 1e-12,
		"opaque volume changed its physical metal thickness");
	grid.SetElectric([](unsigned int n, const Grid::Point &p) { return std::sin(731 * p[0] - 433 * p[1] + 211 * p[2] + n); });
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p) { return std::cos(331 * p[0] + 171 * p[1] - 411 * p[2] + n) / 376.730313668; });
	double expectedE = 0, expectedH = 0;
	for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
	{
		const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
		const double volume = length[0] * length[1] * length[2];
		for (unsigned int n = 0; n < 3; ++n)
		{
			auto a = lo, b = lo;
			for (unsigned int k = 0; k < 3; ++k)
			{
				a[k] += length[k] * (k == n ? 1. / 3 : .5);
				b[k] += length[k] * (k == n ? 2. / 3 : .5);
			}
			const double ha = grid.Magnetic(n, a), hb = grid.Magnetic(n, b);
			if (grid.LeafIsOpaque(leaf))
			{
				Require(ha == 0 && hb == 0 && grid.Electric(n, a) == 0 && grid.Electric(n, b) == 0,
					"excluded conductor exposes spurious bulk fields");
				continue;
			}
			const double h0 = 2 * ha - hb, h1 = 2 * hb - ha;
			expectedH += mu * volume / 4 * (h0 * h0 + h1 * h1);
			for (unsigned int i = 0; i < 2; ++i)
			{
				for (unsigned int j = 0; j < 2; ++j)
				{
					const double e = grid.Evaluate(grid.LeafEdgeVoltage(leaf, n, i, j)) / length[n];
					expectedE += epsilon * volume / 8 * e * e;
				}
			}
		}
	}
	Require(std::abs(grid.ElectricEnergy() / expectedE - 1) < 2e-12 &&
		std::abs(grid.MagneticEnergy() / expectedH - 1) < 2e-12,
		"bulk mass/energy includes excluded conductor volume");
	const auto faces = grid.OpaqueBoundaryFaces();
	Require(faces.size() == 16, "opaque face extraction includes internal faces or misses exposed faces");
	double surfaceArea = 0;
	for (const auto &face : faces)
	{
		const auto length = grid.LeafLengths(face.leaf);
		surfaceArea += length[(face.normal + 1) % 3] * length[(face.normal + 2) % 3];
	}
	Require(std::abs(surfaceArea / (2e-6 + .004 * thickness) - 1) < 2e-12,
		"exposed conductor area does not preserve thickness and sidewalls");
	Grid reference = grid;
	const double dt = grid.StableTimestep() * .8, initial = grid.Energy(dt);
	std::unique_ptr<GPU_Backend> backend;
	std::unique_ptr<GPU_CartesianSubgrid> device;
	if (target)
	{
		backend.reset(GPU_Backend::New(target->backend));
		Require(bool(backend), "opaque-volume device backend unavailable");
		target->Verify(*backend);
		device = backend->CreateCartesianSubgrid(grid);
	}
	double fieldError = 0, peak = 0, balance = 0;
	for (unsigned int step = 0; step < 512; ++step)
	{
		reference.Step(dt, 0);
		if (device)
		{
			device->Step(dt, 0);
			device->ReadFields(grid);
		}
		else
		{
			grid.Step(dt, 0);
		}
		balance = std::max(balance, std::abs(grid.Energy(dt) / initial - 1));
		for (std::size_t leaf = 0; leaf < grid.CellCount(); ++leaf)
		{
			const auto lo = grid.LeafStart(leaf), length = grid.LeafLengths(leaf);
			auto centre = lo;
			for (unsigned int n = 0; n < 3; ++n)
			{
				centre[n] += .5 * length[n];
			}
			for (unsigned int n = 0; n < 3; ++n)
			{
				const double expected = reference.Electric(n, centre);
				peak = std::max(peak, std::abs(expected));
				fieldError = std::max(fieldError, std::abs(grid.Electric(n, centre) - expected));
				const double h = reference.Magnetic(n, centre) * 376.730313668;
				peak = std::max(peak, std::abs(h));
				fieldError = std::max(fieldError, std::abs(grid.Magnetic(n, centre) * 376.730313668 - h));
				if (grid.LeafIsOpaque(leaf))
				{
					Require(grid.Electric(n, centre) == 0 && grid.Magnetic(n, centre) == 0,
						"excluded conductor develops internal fields");
				}
			}
		}
	}
	Require(balance < (target ? 2e-4 : 2e-11) && fieldError / peak < (target ? 2e-4 : 2e-11),
		"opaque-volume fields or passive invariant drifted");
	std::cout << "OPAQUE thickness=" << thickness << " leaves=" << opaqueCount << " boundary_traces=" << surfaceEdges
		<< " exposed_faces=" << faces.size() << " device=" << bool(target) << " field_error=" << fieldError / peak << " balance=" << balance << std::endl;
}
int main(int argc, char **argv)
{
	try
	{
		SubgridDeviceTarget selected;
		const SubgridDeviceTarget *target = nullptr;
		if (argc == 2)
		{
			selected = std::string(argv[1]) == "--reference" ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromFlag(argv[1]);
			target = &selected;
		}
		Control(18e-6, target);
		Control(35e-6, target);
		std::cout << "PASS: excluded physical conductor volume, boundary traces and external-field invariant" << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << std::endl;
		return 1;
	}
}

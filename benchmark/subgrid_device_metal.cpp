/* SPDX-License-Identifier: GPL-3.0-or-later */
// Batch wall time includes completion/readback, never GPU submission time alone.
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_backend.h"
#include "tools/useful.h"
#include <mach-o/dyld.h>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

int main()
{
	try
	{
		SetThreadLimit(1);
		for (std::uint32_t i = 0; i < _dyld_image_count(); ++i)
		{
			const std::string name = _dyld_get_image_name(i);
			if (name.find("libopenEMS") != std::string::npos)
			{
				std::cout << "loaded=" << name << std::endl;
			}
		}
		CartesianSubgrid grid({24, 24, 24}, {.001, .001, .001}, {{6, 6, 6}, {18, 18, 18}},
			8.85418781762e-12, 1.2566370614359173e-6);
		grid.SetElectric([](unsigned int n, const CartesianSubgrid::Point &p)
		{
			return .1 * std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n);
		});
		std::unique_ptr<GPU_Backend> backend(GPU_Backend::New("auto"));
		if (!backend || backend->GetName().find("Metal (") != 0)
		{
			throw std::runtime_error("physical Metal required");
		}
		auto execution = backend->CreateCartesianSubgrid(grid);
		const double dt = .8 * grid.StableTimestep();
		for (unsigned int i = 0; i < 400; ++i)
		{
			execution->Step(dt, 0);
		}
		execution->ReadFields(grid);
		const auto start = std::chrono::steady_clock::now();
		for (unsigned int i = 0; i < 10000; ++i)
		{
			execution->Step(dt, 0);
		}
		execution->ReadFields(grid);
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
		std::cout << "steps=10000 seconds=" << seconds << " bytes=" << execution->StorageBytes()
			<< " energy=" << grid.Energy(dt) << std::endl;
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}

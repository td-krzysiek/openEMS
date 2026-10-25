// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/gpu_backend_domains.h"
#include "FDTD/gpu_plugin.h"
#include "FDTD/operator.h"
#include "openems.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

class EnergySimulation : public openEMS
{
public:
	const Operator& Op() const
	{
		return *FDTD_Op;
	}
};

int main(int argc, char** argv)
{
	try
	{
		if (argc != 2)
		{
			throw std::invalid_argument("usage: test_domain_energy public-fixture.xml");
		}
		EnergySimulation simulation;
		simulation.SetLibraryArguments({ "engine=basic", "numThreads=1" });
		if (!simulation.ParseFDTDSetup(argv[1]) || simulation.SetupFDTD() != 0)
		{
			throw std::runtime_error("fixture setup failed");
		}
		std::vector<std::string> messages;
		const GPU_Plugin* plugin = nullptr;
		for (const auto* candidate : GPU_LoadPlugins(messages))
		{
			if (std::string(candidate->name).find("hip") != std::string::npos &&
			    candidate->DeviceCount() >= 2)
			{
				plugin = candidate;
				break;
			}
		}
		if (!plugin)
		{
			throw std::runtime_error(
			    "this manual test requires a matching HIP plugin and two visible devices");
		}
		const auto& op = simulation.Op();
		unsigned int shape[3];
		for (int n = 0; n < 3; ++n)
		{
			shape[n] = op.GetNumberOfLines(n, true);
		}
		ArrayLib::ArrayNIJK<float> electric, magnetic;
		electric.Init("test_electric", shape);
		magnetic.Init("test_magnetic", shape);
		for (unsigned int c = 0; c < 3; ++c)
		{
			for (unsigned int x = 0; x < shape[0]; ++x)
			{
				for (unsigned int y = 0; y < shape[1]; ++y)
				{
					for (unsigned int z = 0; z < shape[2]; ++z)
					{
						electric(c, x, y, z) = float(1 + c) + .03f * x + .04f * y + .07f * z;
						magnetic(c, x, y, z) = .2f + float(c) + .08f * x + .01f * y + .02f * z;
					}
				}
			}
		}
		for (unsigned int devices : { 1u, 2u })
		{
			std::unique_ptr<GPU_Backend> backend;
			if (devices == 1)
			{
				backend.reset(plugin->NewBackend(0));
			}
			else
			{
				std::vector<std::unique_ptr<GPU_Backend>> domains;
				for (unsigned int n = 0; n < devices; ++n)
				{
					domains.emplace_back(plugin->NewBackend(n));
				}
				backend.reset(new GPU_Backend_Domains(std::move(domains)));
			}
			if (!backend)
			{
				throw std::runtime_error("backend allocation failed");
			}
			backend->Activate();
			if (!backend->Init(&op))
			{
				throw std::runtime_error("backend initialization failed");
			}
			backend->UploadVoltages(electric);
			backend->UploadCurrents(magnetic);
			for (unsigned int nx :
			    { 0u, 1u, shape[0] / 2, shape[0] / 2 + 1, shape[0] - 1, shape[0] })
			{
				unsigned int limit[3] = { nx, shape[1] - 1, shape[2] - 1 };
				double expected_e = 0, expected_h = 0;
				for (unsigned int c = 0; c < 3; ++c)
				{
					for (unsigned int x = 0; x < limit[0]; ++x)
					{
						for (unsigned int y = 0; y < limit[1]; ++y)
						{
							for (unsigned int z = 0; z < limit[2]; ++z)
							{
								const double e = electric(c, x, y, z), h = magnetic(c, x, y, z);
								expected_e += e * e;
								expected_h += h * h;
							}
						}
					}
				}
				double e = 0, h = 0;
				if (!backend->CalcFastEnergy(limit, e, h) ||
				    std::abs(e - expected_e) > 3e-6 * std::max(1., expected_e) ||
				    std::abs(h - expected_h) > 3e-6 * std::max(1., expected_h))
				{
					throw std::runtime_error(
					    "device energy includes a halo or omits/duplicates an owned location");
				}
				std::cout << "PASS devices=" << devices << " owned-prefix=" << nx << " E=" << e
				          << " H=" << h << '\n';
			}
		}
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}

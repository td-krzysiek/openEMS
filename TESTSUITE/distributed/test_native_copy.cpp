// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/gpu_backend_reference.h"
#include "FDTD/gpu_plugin.h"
#include "FDTD/operator.h"
#include "openems.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

class CopySimulation : public openEMS
{
public:
	const Operator& Op() const
	{
		return *FDTD_Op;
	}
};

static void Require(bool condition, const char* message)
{
	if (!condition)
	{
		throw std::runtime_error(message);
	}
}

int main(int argc, char** argv)
{
	try
	{
		Require(argc == 3, "usage: test_domain_copy public-fixture.xml plugin-name|reference");
		CopySimulation simulation;
		simulation.SetLibraryArguments({ "engine=basic", "numThreads=1" });
		Require(simulation.ParseFDTDSetup(argv[1]) && simulation.SetupFDTD() == 0,
		    "fixture setup failed");
		const GPU_Plugin* plugin = nullptr;
		const bool reference = std::string(argv[2]) == "reference";
		if (!reference)
		{
			std::vector<std::string> messages;
			for (const auto* candidate : GPU_LoadPlugins(messages))
			{
				if (std::string(candidate->name).find(argv[2]) != std::string::npos &&
				    candidate->DeviceCount())
				{
					plugin = candidate;
					break;
				}
			}
			Require(plugin, "selected test plugin has no device");
		}
		for (bool domain : { false, true })
		{
			std::unique_ptr<GPU_Backend> backend(
			    reference ? new GPU_Backend_Reference() : plugin->NewBackend(0));
			Require(bool(backend), "backend allocation failed");
			backend->Activate();
			const auto& op = simulation.Op();
			unsigned int shape[3];
			for (unsigned int n = 0; n < 3; ++n)
			{
				shape[n] = op.GetNumberOfLines(n, true);
			}
			Require(domain ? backend->InitDomain(&op, 1, shape[0] - 3) : backend->Init(&op),
			    "initialization failed");
			if (domain)
			{
				shape[0] -= 2;
			}
			const size_t plane = size_t(shape[1]) * shape[2];
			for (bool currents : { false, true })
			{
				ArrayLib::ArrayNIJK<float> expected, actual;
				expected.Init("copy_seed", shape);
				actual.Init("copy_result", shape);
				for (size_t n = 0; n < expected.size(); ++n)
				{
					expected.data()[n] = (currents ? -.25f : .5f) + float(n % 97);
				}
				if (currents)
				{
					backend->UploadCurrents(expected);
				}
				else
				{
					backend->UploadVoltages(expected);
				}
				// Queue dependent copies without an intervening host wait. Cross components
				// and include the exact owned-to-ghost plane used by the native adapter.
				const size_t copies[][3] = { { plane, 0, plane },
					{ 0, expected.size() - plane, plane },
					{ expected.size() - plane, 2 * plane, plane } };
				for (const auto& copy : copies)
				{
					Require(backend->CopyRange(currents, copy[0], copy[1], copy[2]),
					    "valid field copy rejected");
					std::copy(expected.data() + copy[0],
					    expected.data() + copy[0] + copy[2],
					    expected.data() + copy[1]);
				}
				const size_t size = expected.size(), overflow = std::numeric_limits<size_t>::max();
				Require(backend->CopyRange(currents, size, size, 0), "empty end copy rejected");
				Require(backend->CopyRange(currents, 0, 0, size), "identical copy rejected");
				Require(!backend->CopyRange(currents, 0, 1, 2), "overlapping copy accepted");
				Require(!backend->CopyRange(currents, size, 0, 1), "source overflow accepted");
				Require(!backend->CopyRange(currents, 0, size, 1), "destination overflow accepted");
				Require(!backend->CopyRange(currents, 1, 0, overflow), "count overflow accepted");
				Require(!backend->CopyRange(currents, overflow, 0, 1), "offset overflow accepted");
				if (currents)
				{
					backend->DownloadCurrents(actual);
				}
				else
				{
					backend->DownloadVoltages(actual);
				}
				Require(
				    std::equal(expected.data(), expected.data() + expected.size(), actual.data()),
				    "ordered device copy differs or touches another location");
				std::cout << "PASS domain=" << domain << " currents=" << currents << '\n';
			}
		}
	}
	catch (const std::exception& error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}

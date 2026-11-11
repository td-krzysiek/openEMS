// SPDX-License-Identifier: GPL-3.0-or-later
#include "fixture.h"
#include "FDTD/distributed/operator_slice.h"
#include "FDTD/operator.h"
#include "FDTD/engine.h"
#include "openems.h"
#include <iostream>
#include <memory>
#include <fstream>
#include <iomanip>

using namespace openems::distributed;

class Simulation : public openEMS
{
public:
	const Operator& Op() const
	{
		return *FDTD_Op;
	}

	Engine& Fields() const
	{
		return *FDTD_Eng;
	}
};

int main(int argc, char** argv)
{
	if (argc != 2 && argc != 3)
	{
		std::cerr << "Usage: test-operator fixture.xml [coefficient-fixture.txt]\n";
		return 2;
	}
	try
	{
		Simulation simulation;
		simulation.SetLibraryArguments({ "engine=basic", "numThreads=1" });
		if (!simulation.ParseFDTDSetup(argv[1]) || simulation.SetupFDTD() != 0)
		{
			throw std::runtime_error("fixture setup failed");
		}
		const auto& op = simulation.Op();
		auto& engine = simulation.Fields();
		const GridShape shape = {
			op.GetNumberOfLines(0, true), op.GetNumberOfLines(1, true), op.GetNumberOfLines(2, true)
		};
		if (argc == 3)
		{
			PartitionPlan whole(1, shape, { 0, shape.nx });
			const auto& d = whole.Domains()[0];
			auto snapshot = SliceOperator(op, d);
			snapshot.initial_fields = test_fixture::Data(d).initial_fields;
			for (auto& value : snapshot.initial_fields[1])
			{
				value *= .002f;
			}
			std::ofstream file(argv[2]);
			file << "OPENEMS_DOMAIN_OPERATOR_V1\n"
			     << shape.nx << ' ' << shape.ny << ' ' << shape.nz << '\n';
			file << std::setprecision(17) << snapshot.timestep << '\n' << std::setprecision(9);
			for (const auto& a : snapshot.coefficients)
			{
				for (auto value : a)
				{
					file << value << '\n';
				}
			}
			for (const auto& a : snapshot.initial_fields)
			{
				for (auto value : a)
				{
					file << value << '\n';
				}
			}
			if (!file)
			{
				throw std::runtime_error("coefficient fixture export failed");
			}
		}
		PartitionPlan plan(1, shape, { 0, 1, shape.nx / 3, shape.nx - 1, shape.nx });
		std::vector<std::unique_ptr<CPUDomainExecutor>> owned;
		std::vector<DomainExecutor*> executors;
		for (const auto& d : plan.Domains())
		{
			auto data = SliceOperator(op, d);
			data.initial_fields = test_fixture::Data(d).initial_fields;
			for (uint32_t c = 0; c < 3; ++c)
			{
				for (uint32_t x = d.begin; x < d.end; ++x)
				{
					for (uint32_t y = 0; y < shape.ny; ++y)
					{
						for (uint32_t z = 0; z < shape.nz; ++z)
						{
							const size_t i = c * d.OwnedNodes() +
							    size_t(x - d.begin) * shape.Plane() + size_t(y) * shape.nz + z;
							// Scale initial currents to the operator's physical V/I units.
							data.initial_fields[1][i] *= .002f;
							engine.SetVolt(c, x, y, z, data.initial_fields[0][i]);
							engine.SetCurr(c, x, y, z, data.initial_fields[1][i]);
						}
					}
				}
			}
			owned.emplace_back(new CPUDomainExecutor(data));
			executors.push_back(owned.back().get());
		}
		LocalHaloTransport transport;
		LocalCollectiveContext collectives;
		FDTDCoordinator coordinator(plan, transport, collectives, executors);
		RunManager manager(plan, test_fixture::Placements(plan), coordinator);
		manager.Run(120);
		if (!engine.IterateTS(120))
		{
			throw std::runtime_error("basic reference failed");
		}
		double error = 0, peak = 0;
		for (const auto& executor : owned)
		{
			for (int f = 0; f < 2; ++f)
			{
				const auto values =
				    executor->OwnedField(f ? FieldPhase::Current : FieldPhase::Voltage);
				const auto& d = plan.Domains()[executor->Id()];
				for (uint32_t c = 0; c < 3; ++c)
				{
					for (uint32_t x = d.begin; x < d.end; ++x)
					{
						for (uint32_t y = 0; y < shape.ny; ++y)
						{
							for (uint32_t z = 0; z < shape.nz; ++z)
							{
								const size_t i = c * d.OwnedNodes() +
								    size_t(x - d.begin) * shape.Plane() + size_t(y) * shape.nz + z;
								const double reference =
								    f ? engine.GetCurr(c, x, y, z) : engine.GetVolt(c, x, y, z);
								if (!std::isfinite(values[i]) || !std::isfinite(reference))
								{
									throw std::runtime_error("nonfinite operator fields");
								}
								error = std::max(error, std::abs(values[i] - reference));
								peak = std::max(peak, std::abs(reference));
							}
						}
					}
				}
			}
		}
		if (peak < 1e-15 || error / peak > 2e-5)
		{
			throw std::runtime_error("real operator/basic engine mismatch");
		}
		std::cout << "PASS: real global-operator slices, four CPU domains, basic engine error="
		          << error / peak << '\n';
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
}

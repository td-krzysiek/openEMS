// SPDX-License-Identifier: GPL-3.0-or-later
#include "fixture.h"
#include "FDTD/hip/hip_domain_executor.h"
#include <hip/hip_runtime.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <algorithm>

using namespace openems::distributed;
static LocalOperatorData snapshot;
static bool use_snapshot = false;

static void Load(const std::string& path)
{
	std::ifstream in(path);
	std::string magic;
	GridShape shape;
	in >> magic >> shape.nx >> shape.ny >> shape.nz >> snapshot.timestep;
	if (!in || magic != "OPENEMS_DOMAIN_OPERATOR_V1")
	{
		throw std::runtime_error("invalid coefficient fixture");
	}
	shape.Validate();
	if (shape.Nodes() > 2000000)
	{
		throw std::runtime_error("coefficient fixture too large for short test");
	}
	snapshot.domain = { 0, 0, shape.nx, shape };
	for (auto& a : snapshot.coefficients)
	{
		a.resize(3 * shape.Nodes());
		for (auto& v : a)
		{
			in >> v;
		}
	}
	for (auto& a : snapshot.initial_fields)
	{
		a.resize(3 * shape.Nodes());
		for (auto& v : a)
		{
			in >> v;
		}
	}
	if (!in)
	{
		throw std::runtime_error("truncated coefficient fixture");
	}
	in >> std::ws;
	if (!in.eof())
	{
		throw std::runtime_error("trailing coefficient fixture data");
	}
	snapshot.Validate();
	use_snapshot = true;
}

static LocalOperatorData Data(const Domain& d)
{
	if (!use_snapshot)
	{
		return test_fixture::Data(d);
	}
	LocalOperatorData data;
	data.domain = d;
	data.timestep = snapshot.timestep;
	const size_t plane = d.grid.Plane(), global = d.grid.Nodes(), owned = d.OwnedNodes();
	for (size_t n = 0; n < 6; ++n)
	{
		auto& destination = n < 4 ? data.coefficients[n] : data.initial_fields[n - 4];
		const auto& source = n < 4 ? snapshot.coefficients[n] : snapshot.initial_fields[n - 4];
		destination.resize(3 * owned);
		for (uint32_t c = 0; c < 3; ++c)
		{
			std::copy_n(source.data() + c * global + size_t(d.begin) * plane,
			    owned,
			    destination.data() + c * owned);
		}
	}
	data.Validate();
	return data;
}

static double Case(GridShape shape, std::vector<uint32_t> cuts, bool cpu_first)
{
	PartitionPlan p(2, shape, cuts), whole(1, shape, { 0, shape.nx });
	std::vector<std::unique_ptr<HIPDomainExecutor>> gpu;
	std::vector<std::unique_ptr<CPUDomainExecutor>> cpu;
	std::vector<DomainExecutor*> executors;
	for (const auto& d : p.Domains())
	{
		if (cpu_first && d.id == 0)
		{
			cpu.emplace_back(new CPUDomainExecutor(Data(d)));
			executors.push_back(cpu.back().get());
		}
		else
		{
			gpu.emplace_back(new HIPDomainExecutor(Data(d), int(d.id % 2)));
			executors.push_back(gpu.back().get());
		}
	}
	LocalHaloTransport transport;
	LocalCollectiveContext collectives;
	FDTDCoordinator coordinator(p, transport, collectives, executors);
	RunManager manager(p, test_fixture::Placements(p), coordinator);
	manager.Run(120);
	CPUDomainExecutor reference(Data(whole.Domains()[0]));
	LocalHaloTransport rt;
	FDTDCoordinator rc(whole, rt, collectives, { &reference });
	RunManager rm(whole, test_fixture::Placements(whole), rc);
	rm.Run(120);
	std::vector<std::array<std::vector<float>, 2>> fields(p.Domains().size());
	for (const auto& e : gpu)
	{
		fields[e->Id()] = { { e->OwnedField(FieldPhase::Voltage),
			e->OwnedField(FieldPhase::Current) } };
	}
	for (const auto& e : cpu)
	{
		fields[e->Id()] = { { e->OwnedField(FieldPhase::Voltage),
			e->OwnedField(FieldPhase::Current) } };
	}
	const double error = test_fixture::Error(p, reference, fields);
	if (error > 2e-5)
	{
		throw std::runtime_error("HIP/domain six-field comparison failed");
	}
	return error;
}

static void Benchmark(unsigned int domains, unsigned int repetition)
{
	const GridShape shape = { 512, 384, 128 };
	std::vector<uint32_t> cuts = { 0 };
	for (unsigned int n = 1; n <= domains; ++n)
	{
		cuts.push_back(shape.nx * n / domains);
	}
	PartitionPlan p(1, shape, cuts);
	LocalHaloTransport transport;
	LocalCollectiveContext collectives;
	std::vector<std::unique_ptr<HIPDomainExecutor>> gpu;
	std::vector<DomainExecutor*> executors;
	size_t bytes = 0;
	for (const auto& d : p.Domains())
	{
		gpu.emplace_back(new HIPDomainExecutor(Data(d), int(d.id)));
		executors.push_back(gpu.back().get());
		bytes += gpu.back()->AllocatedDeviceBytes();
	}
	FDTDCoordinator coordinator(p, transport, collectives, executors);
	RunManager manager(p, test_fixture::Placements(p), coordinator);
	manager.Start();
	manager.Advance(32);
	const auto start = std::chrono::steady_clock::now();
	manager.Advance(512);
	const double seconds =
	    std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
	manager.Finish();
	std::cout << std::setprecision(12) << "{\"benchmark\":true,\"domains\":" << domains
	          << ",\"repetition\":" << repetition << ",\"nodes\":" << shape.Nodes()
	          << ",\"timed_steps\":512,\"seconds\":" << seconds
	          << ",\"MCells_s\":" << shape.Nodes() * 512. / seconds / 1e6
	          << ",\"device_bytes_total\":" << bytes << "}\n";
}

int main(int argc, char** argv)
{
	try
	{
		int count = 0;
		if (hipGetDeviceCount(&count) != hipSuccess || count != 2)
		{
			throw std::runtime_error("expose exactly two selected GPUs");
		}
		for (int device = 0; device < 2; ++device)
		{
			hipDeviceProp_t prop;
			if (hipGetDeviceProperties(&prop, device) != hipSuccess)
			{
				throw std::runtime_error("cannot query selected GPU");
			}
			char bus[64] = {};
			if (hipDeviceGetPCIBusId(bus, sizeof(bus), device) != hipSuccess)
			{
				throw std::runtime_error("cannot query runtime PCI mapping");
			}
			std::cout << "GPU " << device << ": " << prop.name << " PCI=" << bus << '\n';
		}
		unsigned int cases = 0;
		double largest = 0;
		for (uint32_t nx : { 35, 36 })
		{
			for (uint32_t cut : { 1u, 11u, 18u, nx - 1 })
			{
				largest = std::max(largest, Case({ nx, 19, 13 }, { 0, cut, nx }, false));
				++cases;
			}
		}
		largest = std::max(largest, Case({ 35, 19, 13 }, { 0, 11, 35 }, true));
		++cases;
		std::cout << "{\"correctness\":true,\"cases\":" << cases << ",\"max_error\":" << largest
		          << "}\n";
		if (argc == 3 && std::string(argv[1]) == "--fixture")
		{
			Load(argv[2]);
			const auto shape = snapshot.domain.grid;
			for (const bool cpu : { false, true })
			{
				largest = std::max(largest, Case(shape, { 0, shape.nx / 3, shape.nx }, cpu));
			}
			std::cout << "{\"real_operator\":true,\"max_error\":" << largest << "}\n";
		}
		else if (argc == 2 && std::string(argv[1]) == "--benchmark")
		{
			for (unsigned int rep = 0; rep < 3; ++rep)
			{
				for (unsigned int domains :
				    rep % 2 ? std::vector<unsigned int>{ 2, 1 } : std::vector<unsigned int>{ 1, 2 })
				{
					Benchmark(domains, rep);
				}
			}
		}
		else if (argc != 1)
		{
			throw std::invalid_argument("unknown test argument");
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
}

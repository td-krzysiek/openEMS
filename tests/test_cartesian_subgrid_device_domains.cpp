/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "tools/useful.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <string>

using Grid = CartesianSubgrid;
using Field = CartesianSubgridField;
constexpr double epsilon = 8.85418781762e-12, mu = 1.2566370614359173e-6;

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
	throw std::runtime_error("invalid domain operation was accepted");
}

//! Count explicit host/device traffic; no timing or device fallback machinery.
class CountingBackend : public GPU_Backend
{
  public:
	explicit CountingBackend(GPU_Backend &backend) : real(backend) {}
	void Activate() override { real.Activate(); }
	std::string GetName() const override { return real.GetName(); }
	bool Init(const Operator *op, const GPU_CellModel &model) override { return real.Init(op, model); }
	void *DeviceAlloc(std::size_t bytes, const void *host = nullptr) override
	{
		allocated += bytes;
		return real.DeviceAlloc(bytes, host);
	}
	void DeviceDownload(void *host, const void *dev, std::size_t bytes) override
	{
		downloaded += bytes;
		largestDownload = std::max(largestDownload, bytes);
		real.DeviceDownload(host, dev, bytes);
	}
	void DeviceUpload(void *dev, const void *host, std::size_t bytes) override
	{
		uploaded += bytes;
		real.DeviceUpload(dev, host, bytes);
	}
	float *DeviceField(bool currents) override { return real.DeviceField(currents); }
	void UpdateVoltages(const Engine_Box &box) override { real.UpdateVoltages(box); }
	void UpdateCurrents(const Engine_Box &box) override { real.UpdateCurrents(box); }
	void RunPoints(GPU_PointOp kind, const void *op, uint64_t first, uint64_t stop, bool serial) override
	{
		real.RunPoints(kind, op, first, stop, serial);
	}
	void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT> &field) override { real.Download(currents, field); }
	GPU_Backend *NewSubGridBackend() override { return real.NewSubGridBackend(); }
	void Synchronize() override { real.Synchronize(); }
	std::size_t allocated = 0, downloaded = 0, uploaded = 0, largestDownload = 0;
	GPU_Backend &real;
};

//! Deterministic guard for ambient device state even on a CPU-only test host.
class ContextBackend : public CountingBackend
{
  public:
	ContextBackend(GPU_Backend &backend, std::size_t &active, std::size_t identity)
		: CountingBackend(backend), m_Active(active), m_Identity(identity)
	{
		m_Active = identity;
	}
	void Activate() override
	{
		m_Active = m_Identity;
		real.Activate();
	}
	void *DeviceAlloc(std::size_t bytes, const void *host = nullptr) override
	{
		Check();
		return CountingBackend::DeviceAlloc(bytes, host);
	}
	void DeviceDownload(void *host, const void *device, std::size_t bytes) override
	{
		Check();
		CountingBackend::DeviceDownload(host, device, bytes);
	}
	void DeviceUpload(void *device, const void *host, std::size_t bytes) override
	{
		Check();
		CountingBackend::DeviceUpload(device, host, bytes);
	}
	void RunPoints(GPU_PointOp kind, const void *op, uint64_t first, uint64_t stop, bool serial) override
	{
		Check();
		CountingBackend::RunPoints(kind, op, first, stop, serial);
	}
  private:
	void Check() const
	{
		Require(m_Active == m_Identity, "domain used another backend's ambient context");
	}
	std::size_t &m_Active;
	std::size_t m_Identity;
};

//! A rebuilt external backend may still lack the optional upload capability.
class MissingUploadBackend : public CountingBackend
{
  public:
	using CountingBackend::CountingBackend;
	void DeviceUpload(void *device, const void *host, std::size_t bytes) override
	{
		GPU_Backend::DeviceUpload(device, host, bytes);
	}
};

static Grid Fixture(bool multiple, bool pec)
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
	const auto regions = multiple ?
		std::vector<Grid::Region>{{{{2, 2, 2}, {5, 9, 9}}, {2, 2, 1}},
			{{{7, 3, 3}, {10, 10, 10}}, {1, 2, 4}}} :
		std::vector<Grid::Region>{{{{3, 3, 3}, {9, 9, 9}}, {2, 4, 2}}};
	const auto medium = [&](const Grid::Point &, const Grid::Point &stop)
	{
		return Grid::Medium{epsilon * (stop[0] <= axes[0][6] ? 2.25 : 1), mu};
	};
	const auto metal = [&](unsigned int, const Grid::Point &start, const Grid::Point &stop)
	{
		return pec && start[2] == axes[2][2] && stop[2] == axes[2][2];
	};
	Grid grid(axes, epsilon, mu, regions, false, medium, metal);
	grid.SetElectric([](unsigned int n, const Grid::Point &p)
	{
		return std::sin(1731 * p[0] + 933 * p[1] - 2671 * p[2] + n);
	});
	grid.SetMagnetic([](unsigned int n, const Grid::Point &p)
	{
		return std::cos(933 * p[0] - 1731 * p[1] + 411 * p[2] + n) / 376.730313668;
	});
	grid.SetSources({{grid.EdgeCount() / 2, .1}, {grid.EdgeCount() / 3, -.03}});
	return grid;
}

static double Error(const std::vector<float> &a, const std::vector<float> &b)
{
	Require(a.size() == b.size(), "domain field dimensions differ");
	double delta = 0, scale = 0;
	for (std::size_t i = 0; i < a.size(); ++i)
	{
		Require(std::isfinite(b[i]), "nonfinite domain field");
		delta = std::max(delta, std::abs(double(a[i]) - b[i]));
		scale = std::max(scale, std::abs(double(a[i])));
	}
	return delta / std::max(scale, 1e-300);
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, unsigned int count,
	bool multiple, bool pec, unsigned int steps, bool mixed = false, bool context = false)
{
	auto original = Fixture(multiple, pec), observed = Fixture(multiple, pec), precise = Fixture(multiple, pec);
	const auto &lines = original.BaseLines()[axis];
	std::vector<double> cuts;
	for (unsigned int k = 1; k < count; ++k)
	{
		cuts.push_back(lines.front() + (lines.back() - lines.front()) * k / count);
	}
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(original, axis, cuts);
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<std::unique_ptr<CountingBackend>> counted;
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers;
	std::size_t activeContext = 0;
	for (const auto &domain : plan->Domains())
	{
		target.SelectDomain(backends.size());
		const bool cpuWorker = mixed && backends.size() % 2;
		backends.emplace_back(GPU_Backend::New(cpuWorker ? "reference" : target.backend));
		Require(bool(backends.back()), "selected domain backend unavailable");
		(cpuWorker ? SubgridDeviceTarget{} : target).Verify(*backends.back());
		if (context)
		{
			counted.push_back(std::make_unique<ContextBackend>(*backends.back(), activeContext, backends.size()));
		}
		else
		{
			counted.push_back(std::make_unique<CountingBackend>(*backends.back()));
		}
		workers.push_back(CreateCartesianSubgridWorker(*counted.back(), domain));
	}
	auto transport = CreateCartesianSubgridHostTransport(plan, workers);
	auto *metrics = transport.get();
	auto execution = CreateCartesianSubgridCoordinator(plan, std::move(workers), std::move(transport));
	target.SelectDomain(0);
	std::unique_ptr<GPU_Backend> one(GPU_Backend::New(target.backend));
	Require(bool(one), "selected one-domain backend unavailable");
	target.Verify(*one);
	auto reference = one->CreateCartesianSubgrid(original);
	const double dt = .8 * original.StableTimestep();
	execution->RequireFixedTimestep(dt);
	reference->RequireFixedTimestep(dt);
	Reject([&] { execution->UpdateMagnetic(dt); });
	Reject([&] { execution->Step(2 * dt, 0); });
	Reject([&] { execution->Fields(); });
	for (unsigned int step = 0; step < steps; ++step)
	{
		const double drive = std::exp(-std::pow((double(step) - 20) / 8, 2));
		execution->UpdateElectric(dt, drive);
		if (!step)
		{
			Reject([&] { execution->ReadFields(observed); });
			Reject([&] { execution->UpdateElectric(dt, drive); });
			Reject([&] { execution->UpdateMagnetic(dt * .5); });
		}
		execution->UpdateMagnetic(dt);
		reference->Step(dt, drive);
		precise.Step(dt, drive);
	}
	std::size_t traffic = 0, largestRead = 0;
	for (const auto &backend : counted)
	{
		traffic += backend->downloaded + backend->uploaded;
		largestRead = std::max(largestRead, backend->largestDownload);
	}
	Require(traffic == metrics->HostTrafficBytes(), "unaccounted per-step device transfer");
	const auto fullBytes = std::max(plan->ElectricCount(), plan->MagneticCount()) * sizeof(float);
	Require(largestRead < fullBytes, "halo exchange downloaded a full global field");
	if (count > 1)
	{
		Require(traffic > 0 && largestRead > 0, "no actual device halo exchange");
	}
	else
	{
		Require(traffic == 0, "one domain exchanged unnecessary halos");
	}
	execution->ReadFields(observed);
	reference->ReadFields(original);
	const auto a = original.DeviceData(), b = observed.DeviceData(), c = precise.DeviceData();
	const double singleError = std::max(Error(a.electric, b.electric), Error(a.magnetic, b.magnetic));
	const double doubleError = std::max(Error(c.electric, b.electric), Error(c.magnetic, b.magnetic));
	Require(singleError < 2e-6, "partition changed one-domain FP32 fields");
	Require(doubleError < 2e-4, "partition/double field error exceeded 2e-4");
	Require(std::abs(observed.Energy(dt) / original.Energy(dt) - 1) < 2e-6,
		"owned-state energy reconstruction differs");
	std::cout << "axis=" << axis << " domains=" << count << " multiple=" << multiple << " PEC=" << pec
		<< " mixed=" << mixed << " context=" << context << " steps=" << steps << " single_error=" << singleError << " double_error=" << doubleError
		<< " halo_bytes=" << traffic << " largest_download=" << largestRead
		<< " storage=" << execution->StorageBytes() << std::endl;
}

//! A lost completion must fail the run, not advance with stale E ghosts.
class MissingElectricExchange : public CartesianSubgridDomainTransport
{
  public:
	explicit MissingElectricExchange(std::unique_ptr<CartesianSubgridDomainTransport> transport)
		: m_Transport(std::move(transport)) {}
	void Exchange(Field field, std::size_t epoch) override
	{
		if (field != Field::Electric)
		{
			m_Transport->Exchange(field, epoch);
		}
	}
	std::uint64_t HostTrafficBytes() const override { return m_Transport->HostTrafficBytes(); }
  private:
	std::unique_ptr<CartesianSubgridDomainTransport> m_Transport;
};

static void LostExchange(const SubgridDeviceTarget &target)
{
	auto grid = Fixture(false, false);
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(grid, 0,
		std::vector<double>{.5 * grid.BaseLines()[0].back()});
	std::vector<std::unique_ptr<GPU_Backend>> backends;
	std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers;
	for (const auto &domain : plan->Domains())
	{
		target.SelectDomain(backends.size());
		backends.emplace_back(GPU_Backend::New(target.backend));
		Require(bool(backends.back()), "selected lost-exchange backend unavailable");
		target.Verify(*backends.back());
		workers.push_back(CreateCartesianSubgridWorker(*backends.back(), domain));
	}
	auto broken = std::make_unique<MissingElectricExchange>(CreateCartesianSubgridHostTransport(plan, workers));
	auto execution = CreateCartesianSubgridCoordinator(plan, std::move(workers), std::move(broken));
	const double dt = .8 * grid.StableTimestep();
	Reject([&] { execution->Step(dt, 0); });
	Reject([&] { execution->Step(dt, 0); });
	Reject([&] { execution->ReadFields(grid); });
	std::cout << "PASS: lost E exchange and poisoned-coordinator continuation rejected" << std::endl;
}

static void Guards(const SubgridDeviceTarget &target)
{
	auto grid = Fixture(false, false);
	const double mid = .5 * grid.BaseLines()[0].back();
	Reject([&] { CartesianSubgridDomainPlan invalid(grid, 3, {mid}); });
	Reject([&] { CartesianSubgridDomainPlan invalid(grid, 0, {mid, mid}); });
	Reject([&] { CartesianSubgridDomainPlan invalid(grid, 0, {std::numeric_limits<double>::quiet_NaN()}); });
	auto plan = std::make_shared<CartesianSubgridDomainPlan>(grid, 0, std::vector<double>{mid});
	std::unique_ptr<GPU_Backend> backend(GPU_Backend::New(target.backend));
	Require(bool(backend), "selected guard backend unavailable");
	target.Verify(*backend);
	auto data = plan->Domains().front();
	data.fields.sourceEdge.push_back(static_cast<std::uint32_t>(data.ownedElectric));
	data.fields.sourceWeight.push_back(1);
	data.fields.sourceOffStart.push_back(0);
	Reject([&] { CreateCartesianSubgridDomainExecution(*backend, data); });
	const CartesianSubgridDomainPlan expanded(grid, 0, {mid},
		CartesianSubgridDomainPlan::CurlPreparation::ExpandedOracle);
	data = expanded.Domains().front();
	data.fields.transposeColumn.front() = static_cast<std::uint32_t>(data.fields.magnetic.size());
	Reject([&] { CreateCartesianSubgridDomainExecution(*backend, data); });
	CountingBackend allocationGuard(*backend);
	const auto malformed = [&](const std::function<void(CartesianSubgridDomainData &)> &change)
	{
		auto invalid = plan->Domains().front();
		Require(invalid.operators.has_value(), "missing compressed negative-control fixture");
		change(invalid);
		const auto before = allocationGuard.allocated;
		Reject([&] { CreateCartesianSubgridDomainExecution(allocationGuard, invalid); });
		Require(allocationGuard.allocated == before, "malformed compressed domain allocated GPU memory");
	};
	malformed([](auto &d) { d.operators.reset(); });
	malformed([](auto &d) { d.fields.streamCurls = false; });
	malformed([](auto &d) { d.fields.transposeColumn.push_back(0); });
	malformed([](auto &d) { d.operators->geometry.push_back(0); });
	malformed([](auto &d) { d.operators->transpose.base.pop_back(); });
	malformed([](auto &d) { d.operators->transpose.weight.clear(); });
	malformed([](auto &d) { d.operators->curl.start.front() = 1; });
	malformed([](auto &d) { d.operators->curl.start.back() = 0; });
	malformed([](auto &d) { d.operators->curl.pattern.front() = std::numeric_limits<std::uint32_t>::max(); });
	malformed([](auto &d) { d.operators->curl.offset.front() = std::numeric_limits<std::uint32_t>::max(); });
	malformed([](auto &d) { d.operators->curl.weight.front() = std::numeric_limits<float>::quiet_NaN(); });
	malformed([](auto &d) { d.operators->curl.residual.front() = std::numeric_limits<float>::denorm_min(); });
	malformed([](auto &d)
	{
		auto &curl = d.operators->curl;
		for (std::size_t row = 0; row < curl.base.size(); ++row)
		{
			const auto pattern = curl.pattern[row];
			if (curl.start[pattern] < curl.start[pattern + 1])
			{
				curl.base[row] = std::numeric_limits<std::uint32_t>::max();
				return;
			}
		}
		throw std::runtime_error("no nonempty compressed row for overflow control");
	});
	std::cout << "PASS: 13 malformed compressed domains rejected before GPU allocation" << std::endl;
	data = plan->Domains().front();
	Require(!data.fields.lower.empty(), "missing factor negative-control fixture");
	data.fields.lower.front() = -1;
	Reject([&] { CreateCartesianSubgridDomainExecution(*backend, data); });
	data = plan->Domains().front();
	data.fields.ledger[0] = 1;
	Reject([&] { CreateCartesianSubgridDomainExecution(*backend, data); });
	const auto halo = std::find_if(plan->Domains().begin(), plan->Domains().end(),
		[](const CartesianSubgridDomainData &d) { return d.ownedMagnetic < d.fields.magnetic.size(); });
	Require(halo != plan->Domains().end(), "missing magnetic ghost guard fixture");
	auto raw = CreateCartesianSubgridDomainExecution(*backend, *halo);
	Reject([&] { raw->UpdateElectric(.8 * grid.StableTimestep(), 0); });
	Reject([&] { raw->AcknowledgeGhosts(true, 0); });
	Reject([&] { raw->AcknowledgeGhosts(false, 1); });
	Reject([&] { raw->ReadFields(grid); });
	auto worker = CreateCartesianSubgridWorker(*backend, *halo);
	std::vector<float> nonfinite(worker->Ghosts(Field::Magnetic), std::numeric_limits<float>::quiet_NaN());
	Reject([&] { worker->InstallGhosts(Field::Magnetic, 0, nonfinite.data(), nonfinite.size()); });
	Reject([&] { worker->Electric(.8 * grid.StableTimestep(), 0, 1); });
	MissingUploadBackend unavailable(*backend);
	Reject([&] { CreateCartesianSubgridWorker(unavailable, *halo); });
	Require(unavailable.allocated == 0, "unsupported domain backend allocated fields before rejection");
	// The optional operation must not disable the unchanged one-domain path.
	auto compatible = Fixture(false, false);
	auto legacy = unavailable.CreateCartesianSubgrid(compatible);
	legacy->Step(.8 * compatible.StableTimestep(), 0);
	legacy->ReadFields(compatible);
	auto unownedHistory = grid;
	unownedHistory.ImportDeviceLedger(1, 0);
	Reject([&] { CartesianSubgridDomainPlan invalid(unownedHistory, 0, {mid}); });
	const auto &axes = grid.BaseLines();
	grid.SetSources({});
	grid.SetResistors({{{axes[0][1], axes[1][1], axes[2][1]},
		{axes[0][1], axes[1][1], axes[2][2]}, 50, 0}}, .7 * grid.StableTimestep());
	const CartesianSubgridDomainPlan ports(grid, 0, {mid});
	Require(ports.PortTimestep() == .7 * grid.StableTimestep(), "port recurrence step was lost");
	std::cout << "PASS: domain cuts and malformed local uploads reject; coupled port plan retains timestep" << std::endl;
}

int main(int argc, char **argv)
{
	try
	{
		SetThreadLimit(1);
		Require(argc <= 2, "usage: test_cartesian_subgrid_device_domains [--metal|--cuda|--cuda-two|--cuda-two-reverse]");
		const SubgridDeviceTarget target = argc == 2 ? SubgridDeviceTarget::FromDomainFlag(argv[1]) : SubgridDeviceTarget{};
		Control(target, 0, 1, false, false, 1500);
		for (unsigned int axis = 0; axis < 3; ++axis)
		{
			Control(target, axis, 2 + axis, false, false, 1500);
			Control(target, axis, 2 + axis, true, axis == 2, 1500);
		}
		if (target.backend != "reference")
		{
			Control(target, 0, 2, false, false, 1500, true);
		}
		Control(target, 2, 3, true, true, 128, false, true);
		Guards(target);
		LostExchange(target);
		std::cout << "PASS: coordinator/local domains and boundary-only transport" << std::endl;
	}
	catch (const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
	return 0;
}

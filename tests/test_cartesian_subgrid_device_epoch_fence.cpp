/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/cartesian_subgrid_domains.h"
#include "FDTD/gpu_backend.h"
#include "subgrid_device_target.h"
#include "subgrid_fence_backend.h"
#include "tools/useful.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>
using Grid = CartesianSubgrid;
static void Require(bool condition, const char *message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}
template<class Action> static void Reject(Action action)
{
    try
    {
        action();
    }
    catch (const std::logic_error &)
    {
        return;
    }
    throw std::runtime_error("invalid epoch fence accepted");
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

class FaultWorker : public CartesianSubgridDomainWorker
{
  public:
    FaultWorker(const CartesianSubgridDomainData &data, unsigned int fault, unsigned int &calls)
        : m_Data(data), m_Fault(fault), m_Calls(calls) {}
    void Electric(double, double, std::size_t) override { throw std::logic_error("unused fake electric"); }
    void Magnetic(double, std::size_t) override { throw std::logic_error("unused fake magnetic"); }
    void RequireFixedTimestep(double) override {}
    std::size_t Owned(CartesianSubgridField field) const override
    {
        return field == CartesianSubgridField::Electric ? m_Data.ownedElectric : m_Data.ownedMagnetic;
    }
    std::size_t Exports(CartesianSubgridField field) const override
    {
        return field == CartesianSubgridField::Electric ? m_Data.sendElectric.size() : m_Data.sendMagnetic.size();
    }
    std::size_t Ghosts(CartesianSubgridField field) const override
    {
        return (field == CartesianSubgridField::Electric ? m_Data.globalElectric.size() : m_Data.globalMagnetic.size()) - Owned(field);
    }
    void ReadExports(CartesianSubgridField, std::size_t, float *, std::size_t) override { throw std::logic_error("unused fake export"); }
    void InstallGhosts(CartesianSubgridField, std::size_t, const float *, std::size_t) override { throw std::logic_error("unused fake upload"); }
    void ReadOwned(std::size_t, float *, std::size_t, float *, std::size_t) override { throw std::logic_error("unused fake fields"); }
    std::array<double, 2> ReadLedger(std::size_t) override { return {}; }
    std::size_t StorageBytes() const override { return 0; }
    std::size_t FenceObservationEpoch(std::size_t epoch) override
    {
        ++m_Calls;
        if (m_Fault == 1)
        {
            return epoch + 1;
        }
        if (m_Fault == 2)
        {
            throw std::logic_error("injected fence failure");
        }
        if (m_Fault == 3)
        {
            return CartesianSubgridDomainWorker::FenceObservationEpoch(epoch);
        }
        return epoch;
    }

  private:
    const CartesianSubgridDomainData &m_Data;
    unsigned int m_Fault;
    unsigned int &m_Calls;
};

static void Faults()
{
    auto grid = Fixture(false);
    const auto &lines = grid.BaseLines()[0];
    for (unsigned int fault = 0; fault < 4; ++fault)
    {
        auto plan = std::make_shared<CartesianSubgridDomainPlan>(grid, 0, std::vector<double>{lines[3], lines[5]});
        std::array<unsigned int, 3> calls{};
        std::vector<std::unique_ptr<CartesianSubgridDomainWorker>> workers;
        for (unsigned int k = 0; k < 3; ++k)
        {
            workers.push_back(std::make_unique<FaultWorker>(plan->Domains()[k], k == 1 ? fault : 0, calls[k]));
        }
        auto transport = CreateCartesianSubgridHostTransport(plan, workers);
        auto execution = CreateCartesianSubgridCoordinator(plan, std::move(workers), std::move(transport));
        if (!fault)
        {
            Require(execution->FenceObservationEpoch(0) == 3 && calls == std::array<unsigned int, 3>{1, 1, 1},
                "valid acknowledgements were not collected");
        }
        else
        {
            Reject([&] { execution->FenceObservationEpoch(0); });
            Require(calls == std::array<unsigned int, 3>{1, 1, 0}, "mixed/failed fence published partial acknowledgement");
            Reject([&] { execution->FenceObservationEpoch(0); });
            Reject([&] { execution->Step(.7 * grid.StableTimestep(), 0); });
        }
    }
    std::cout << "EPOCH_FENCE_FAULTS wrong_epoch=1 worker_failure=1 unsupported=1 poisoned=1\n";
}

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool ports)
{
    auto original = Fixture(ports), singleGrid = original, splitGrid = original, controlGrid = original;
    const double dt = .7 * original.StableTimestep();
    target.SelectDomain(0);
    auto backend = std::make_unique<SubgridFenceBackend>(
        std::unique_ptr<GPU_Backend>(GPU_Backend::New(target.backend)));
    Require(bool(backend), "fence backend missing");
    target.Verify(*backend);
    auto single = backend->CreateCartesianSubgrid(singleGrid), control = backend->CreateCartesianSubgrid(controlGrid);
    std::vector<std::unique_ptr<SubgridFenceBackend>> backends;
    std::vector<GPU_Backend *> pointers;
    for (unsigned int k = 0; k < 3; ++k)
    {
        target.SelectDomain(k);
        backends.push_back(std::make_unique<SubgridFenceBackend>(
            std::unique_ptr<GPU_Backend>(GPU_Backend::New(target.backend))));
        Require(bool(backends.back()), "fence domain backend missing");
        target.Verify(*backends.back());
        pointers.push_back(backends.back().get());
    }
    const auto &lines = original.BaseLines()[axis];
    auto plan = std::make_shared<CartesianSubgridDomainPlan>(splitGrid, axis, std::vector<double>{lines[3], lines[5]});
    auto split = CreateCartesianSubgridDomains(plan, pointers);
    for (auto *execution : {single.get(), split.get(), control.get()})
    {
        execution->AddObservationDFT(true, {0, static_cast<std::uint32_t>(original.EdgeCount() / 2)}, 1, 4096);
        execution->AddObservationDFT(false, {0, static_cast<std::uint32_t>(original.FaceCount() / 2)}, 1, 4096);
    }
    const auto singleHost = single->HostObservationStorageBytes(), splitHost = split->HostObservationStorageBytes();
    auto fence = [&](std::size_t epoch)
    {
        const auto singleBefore = backend->completedWaits;
        std::array<std::size_t, 3> before{};
        for (std::size_t k = 0; k < backends.size(); ++k)
        {
            before[k] = backends[k]->completedWaits;
        }
        Require(single->FenceObservationEpoch(epoch) == 1 && split->FenceObservationEpoch(epoch) == 3,
            "complete fence count differs");
        Require(backend->completedWaits == singleBefore + 1,
            "single-device acknowledgement skipped its backend completion wait");
        for (std::size_t k = 0; k < backends.size(); ++k)
        {
            Require(backends[k]->completedWaits == before[k] + 1,
                "domain acknowledgement skipped its backend completion wait");
        }
    };
    fence(0);
    for (unsigned int step = 0; step < 64; ++step)
    {
        const double amplitude = .03 * std::sin(.17 * step);
        for (auto *execution : {single.get(), split.get()})
        {
            execution->UpdateElectric(dt, amplitude);
            Reject([&] { execution->FenceObservationEpoch(step); });
            execution->UpdateMagnetic(dt);
        }
        control->Step(dt, amplitude);
        if (step % 8 == 7)
        {
            for (auto *execution : {single.get(), split.get(), control.get()})
            {
                execution->AccumulateObservationDFT(0, {{1e-14f, -.3e-14f}});
                execution->AccumulateObservationDFT(1, {{.2e-14f, 1e-14f}});
            }
            Reject([&] { single->FenceObservationEpoch(step); });
            Reject([&] { split->FenceObservationEpoch(step); });
            for (unsigned int repeat = 0; repeat < 2; ++repeat)
            {
                fence(step + 1);
            }
            Require(single->ObservationDFTDownloadBytes() == 0 && split->ObservationDFTDownloadBytes() == 0 &&
                single->HostObservationStorageBytes() == singleHost && split->HostObservationStorageBytes() == splitHost,
                "fence downloaded/materialized state or DFT");
        }
    }
    single->ReadFields(singleGrid);
    split->ReadFields(splitGrid);
    control->ReadFields(controlGrid);
    const auto a = singleGrid.DeviceStateData(), b = splitGrid.DeviceStateData(), c = controlGrid.DeviceStateData();
    Require(a.electric == c.electric && a.magnetic == c.magnetic, "fence changed single-device fields");
    double delta = 0, scale = 0;
    for (std::size_t k = 0; k < a.electric.size(); ++k)
    {
        delta = std::max(delta, double(std::abs(a.electric[k] - b.electric[k])));
        scale = std::max(scale, double(std::abs(a.electric[k])));
    }
    Require(delta < 2e-5 * scale, "fence changed domain fields");
    double dftError = 0;
    for (unsigned int field = 0; field < 2; ++field)
    {
        std::vector<std::complex<double>> x, y, z;
        single->ReadObservationDFT(field, 0, x);
        split->ReadObservationDFT(field, 0, y);
        control->ReadObservationDFT(field, 0, z);
        Require(x == z, "fence reset or changed single-device DFT");
        for (std::size_t k = 0; k < x.size(); ++k)
        {
            dftError = std::max(dftError, std::abs(x[k] - y[k]));
        }
        double norm = 0;
        for (auto value : x)
        {
            norm = std::max(norm, std::abs(value));
        }
        Require(dftError < 2e-5 * norm, "fence changed domain DFT");
        dftError = 0;
    }
    std::cout << "EPOCH_FENCE axis=" << axis << " ports=" << ports
        << " fences=17 single_workers=1 domain_workers=3 exact_single_fields=1 exact_single_DFT=1 no_fence_download=1 backend_completion_wait_verified=1\n";
}

int main(int argc, char **argv)
{
    try
    {
        SetThreadLimit(1);
        Require(argc <= 2, "unexpected fence arguments");
        const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
        Faults();
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
            Control(target, axis, false);
            Control(target, axis, true);
        }
        std::cout << "PASS: collective completed-epoch fences, failure guards and unchanged resumed fields/DFT\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

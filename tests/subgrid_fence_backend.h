/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef TEST_SUBGRID_FENCE_BACKEND_H
#define TEST_SUBGRID_FENCE_BACKEND_H

#include "FDTD/gpu_backend.h"
#include <memory>
#include <stdexcept>

// Keep the base CreateCartesianSubgrid factory: its execution must retain this
// wrapper, so a later blocking field read cannot hide a missing fence wait.
class SubgridFenceBackend : public GPU_Backend
{
public:
    explicit SubgridFenceBackend(std::unique_ptr<GPU_Backend> backend)
        : m_Backend(std::move(backend))
    {
        if (!m_Backend)
        {
            throw std::invalid_argument("missing fence test backend");
        }
    }
    void Activate() override
    {
        m_Backend->Activate();
    }
    std::string GetName() const override
    {
        return m_Backend->GetName();
    }
    bool Init(const Operator *op, const GPU_CellModel &model) override
    {
        return m_Backend->Init(op, model);
    }
    void *DeviceAlloc(size_t bytes, const void *host = nullptr) override
    {
        return m_Backend->DeviceAlloc(bytes, host);
    }
    void DeviceDownload(void *host, const void *device, size_t bytes) override
    {
        m_Backend->DeviceDownload(host, device, bytes);
    }
    void DeviceUpload(void *device, const void *host, size_t bytes) override
    {
        m_Backend->DeviceUpload(device, host, bytes);
    }
    float *DeviceField(bool currents) override
    {
        return m_Backend->DeviceField(currents);
    }
    void UpdateVoltages(const Engine_Box &box) override
    {
        m_Backend->UpdateVoltages(box);
    }
    void UpdateCurrents(const Engine_Box &box) override
    {
        m_Backend->UpdateCurrents(box);
    }
    void RunPoints(GPU_PointOp kind, const void *op, uint64_t begin,
                   uint64_t end, bool serial) override
    {
        m_Backend->RunPoints(kind, op, begin, end, serial);
    }
    void Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT> &field) override
    {
        m_Backend->Download(currents, field);
    }
    void Synchronize() override
    {
        m_Backend->Synchronize();
        ++completedWaits;
    }
    GPU_Backend *NewSubGridBackend() override
    {
        return new SubgridFenceBackend(
            std::unique_ptr<GPU_Backend>(m_Backend->NewSubGridBackend()));
    }
    size_t completedWaits = 0;

private:
    std::unique_ptr<GPU_Backend> m_Backend;
};
#endif

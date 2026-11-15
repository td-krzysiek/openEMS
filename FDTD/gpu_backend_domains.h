// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_GPU_BACKEND_DOMAINS_H
#define OPENEMS_GPU_BACKEND_DOMAINS_H
#include "gpu_backend.h"
#include <memory>

// Native Engine adapter; partition/exchange policy remains in distributed/.
class GPU_Backend_Domains : public GPU_Backend
{
public:
	explicit GPU_Backend_Domains(
	    std::vector<std::unique_ptr<GPU_Backend>> devices, std::vector<uint32_t> cuts = {});
	~GPU_Backend_Domains() override;
	std::string GetName() const override;
	bool Init(const Operator*, const GPU_CellModel&) override;
	void UpdateVoltages(const Engine_Box&) override;
	void UpdateCurrents(const Engine_Box&) override;
	void FinishVoltagePhase() override;
	void FinishCurrentPhase() override;
	bool HasDomainPhases() const override
	{
		return true;
	}
	void Download(bool, ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	void* DeviceAlloc(size_t, const void*) override;
	void DeviceDownload(void*, const void*, size_t) override;
	float* DeviceField(bool) override;
	void RunPoints(GPU_PointOp, const void*, uint64_t, uint64_t, bool) override;
	bool UploadRange(bool, size_t, size_t, const FDTD_FLOAT*) override;
	bool DownloadRange(bool, size_t, size_t, FDTD_FLOAT*) override;
	bool CalcFastEnergy(
	    const unsigned int numNodes[3], double& electric, double& magnetic) override;
	void Synchronize() override;
	GPU_Extension* DomainExtension(Engine_Extension*, const unsigned int[3]) override;
	GPU_Backend* NewSubGridBackend() override;

private:
	struct Impl;
	std::unique_ptr<Impl> d;
	void Upload(bool, const ArrayLib::ArrayNIJK<FDTD_FLOAT>&);
};
#endif

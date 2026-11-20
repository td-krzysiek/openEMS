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
	void* DeviceAlloc(size_t, const void*) override {throw std::logic_error("allocate on a child domain");}
	void DeviceDownload(void*, const void*, size_t) override {throw std::logic_error("use domain range download");}
	float* DeviceField(bool) override {throw std::logic_error("domain fields are disjoint");}
	void RunPoints(GPU_PointOp, const void*, uint64_t, uint64_t, bool, const CM_PointFilter&) override {throw std::logic_error("dispatch to a child domain");}
	void Upload(bool, const ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	void Download(bool, ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	void UpdateVoltages() override;
	void UpdateCurrents() override;
	void FinishVoltagePhase() override;
	void FinishCurrentPhase() override;
	bool HasDomainPhases() const override
	{
		return true;
	}
	bool DownloadRange(bool, size_t, size_t, FDTD_FLOAT*) override;
	bool CalcFastEnergy(
	    const unsigned int numNodes[3], double& electric, double& magnetic) override;
	void Synchronize() override;
	GPU_Extension* CreateExtension(Engine_Extension*, Engine*) override;
	GPU_Backend* NewSubGridBackend() override;

private:
	struct Impl;
	std::unique_ptr<Impl> d;
};
#endif

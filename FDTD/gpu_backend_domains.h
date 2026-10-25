// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_GPU_BACKEND_DOMAINS_H
#define OPENEMS_GPU_BACKEND_DOMAINS_H
#include "gpu_backend.h"
#include <memory>

// Native Engine adapter; partition/exchange policy remains in distributed/.
class GPU_Backend_Domains : public GPU_Backend
{
public:
	explicit GPU_Backend_Domains(std::vector<std::unique_ptr<GPU_Backend>> devices);
	~GPU_Backend_Domains() override;
	std::string GetName() const override;
	bool Init(const Operator*) override;
	void UpdateVoltages() override;
	void UpdateCurrents() override;
	void FinishVoltagePhase() override;
	void FinishCurrentPhase() override;
	bool HasDomainPhases() const override
	{
		return true;
	}
	void DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	void DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	void UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	void UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>&) override;
	bool DownloadRange(bool, size_t, size_t, FDTD_FLOAT*) override;
	bool CalcFastEnergy(
	    const unsigned int numNodes[3], double& electric, double& magnetic) override;
	void Synchronize() override;
	GPU_Extension* CreateExtension(Engine_Extension*, Engine*) override;
	GPU_Backend* NewSubGridBackend() override;

private:
	struct Impl;
	std::unique_ptr<Impl> d;
	void Download(bool, ArrayLib::ArrayNIJK<FDTD_FLOAT>&);
	void Upload(bool, const ArrayLib::ArrayNIJK<FDTD_FLOAT>&);
};
#endif

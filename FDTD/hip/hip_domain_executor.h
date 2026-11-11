// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_HIP_DOMAIN_EXECUTOR_H
#define OPENEMS_HIP_DOMAIN_EXECUTOR_H
#include "FDTD/distributed/cpu_executor.h"

namespace openems
{
namespace distributed
{

class HIPDomainExecutor : public DomainExecutor
{
public:
	// Native device ordinal belongs to this adapter only; resolve UUIDs at setup.
	HIPDomainExecutor(const LocalOperatorData&, int device);
	~HIPDomainExecutor();
	uint32_t Id() const override;
	double Timestep() const override;
	CompletionPtr Update(FieldPhase, uint64_t) override;
	CompletionPtr Pack(const HaloRecipe&, HostBuffer) override;
	CompletionPtr Unpack(const HaloRecipe&, HostBuffer) override;
	std::vector<float> OwnedField(FieldPhase) const;
	std::string DeviceName() const;
	size_t AllocatedDeviceBytes() const;

private:
	struct Impl;
	std::shared_ptr<Impl> m_Impl;
};

}
}
#endif

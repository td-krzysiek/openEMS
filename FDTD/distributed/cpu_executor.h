// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_CPU_EXECUTOR_H
#define OPENEMS_DISTRIBUTED_CPU_EXECUTOR_H
#include "runtime.h"
#include <array>

namespace openems
{
namespace distributed
{

// Final global-operator coefficients sliced by ownership, component-major.
// c[0..3] are vv,vi,ii,iv; initial fields contain only owned samples.
struct LocalOperatorData
{
	Domain domain;
	double timestep;
	std::array<std::vector<float>, 4> coefficients;
	std::array<std::vector<float>, 2> initial_fields;
	std::vector<std::string> required_extensions;
	void Validate() const;
};

class CPUDomainExecutor : public DomainExecutor
{
public:
	explicit CPUDomainExecutor(const LocalOperatorData&);

	uint32_t Id() const override
	{
		return m_Data.domain.id;
	}

	double Timestep() const override
	{
		return m_Data.timestep;
	}

	CompletionPtr Update(FieldPhase, uint64_t) override;
	CompletionPtr Pack(const HaloRecipe&, HostBuffer) override;
	CompletionPtr Unpack(const HaloRecipe&, HostBuffer) override;
	std::vector<float> OwnedField(FieldPhase) const;

private:
	LocalOperatorData m_Data;
	std::vector<float> m_Voltage, m_Current;
	void Voltage();
	void Current();
};

}
}
#endif

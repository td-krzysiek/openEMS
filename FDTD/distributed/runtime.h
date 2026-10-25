// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_RUNTIME_H
#define OPENEMS_DISTRIBUTED_RUNTIME_H

#include "domain_plan.h"
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace openems
{
namespace distributed
{

class Completion
{
public:
	virtual ~Completion()
	{
	}

	virtual bool Ready() = 0; // success or throws on failure
	virtual void Wait() = 0;
};

typedef std::shared_ptr<Completion> CompletionPtr;
CompletionPtr Completed();

// First wire contract: packed tangential FP32 host samples, component-major.
// Adapters own native device storage and must bridge their own queue events.
typedef std::shared_ptr<std::vector<float>> HostBuffer;

class HaloTransport
{
public:
	virtual ~HaloTransport()
	{
	}

	virtual CompletionPtr Receive(const MessageId& id, HostBuffer destination) = 0;
	virtual CompletionPtr Send(const MessageId& id, HostBuffer source) = 0;
	virtual void Progress() = 0;
	virtual void Fail(const std::string& reason) = 0;
};

class LocalHaloTransport : public HaloTransport
{
public:
	explicit LocalHaloTransport(unsigned int delayed_polls = 0);
	~LocalHaloTransport();
	CompletionPtr Receive(const MessageId&, HostBuffer) override;
	CompletionPtr Send(const MessageId&, HostBuffer) override;
	void Progress() override;
	void Fail(const std::string&) override;

private:
	struct Impl;
	std::shared_ptr<Impl> m_Impl;
};

class CollectiveContext
{
public:
	virtual ~CollectiveContext()
	{
	}

	virtual double Minimum(double value) = 0;
	virtual double Sum(double value) = 0;
	virtual bool All(bool value) = 0;
	virtual void Fail(const std::string& reason) = 0;
};

class LocalCollectiveContext : public CollectiveContext
{
public:
	double Minimum(double) override;
	double Sum(double) override;

	bool All(bool value) override
	{
		return value;
	}

	void Fail(const std::string&) override;
};

class DomainExecutor
{
public:
	virtual ~DomainExecutor()
	{
	}

	virtual uint32_t Id() const = 0;
	virtual double Timestep() const = 0;
	virtual CompletionPtr Update(FieldPhase phase, uint64_t step) = 0;
	virtual CompletionPtr Pack(const HaloRecipe&, HostBuffer) = 0;
	virtual CompletionPtr Unpack(const HaloRecipe&, HostBuffer) = 0;

	virtual void Observe(FieldPhase phase, uint64_t step)
	{
		(void)phase;
		(void)step;
	}
};

class ExchangeResult
{
public:
	ExchangeResult(const ExchangeResult&) = default;
	ExchangeResult& operator=(const ExchangeResult&) = default;

	const ExchangeEpoch& Epoch() const
	{
		return m_Epoch;
	}

	const std::vector<CompletionPtr>& SendRelease() const
	{
		return m_SendRelease;
	}

	const std::vector<CompletionPtr>& ReceiveReady() const
	{
		return m_ReceiveReady;
	}

private:
	friend class HaloExchange;

	ExchangeResult() : m_Epoch{ 0, 0, FieldPhase::Current }
	{
	}

	ExchangeEpoch m_Epoch;
	std::vector<CompletionPtr> m_SendRelease, m_ReceiveReady;
};

class HaloExchange
{
public:
	HaloExchange(const PartitionPlan&, HaloTransport&, const std::vector<DomainExecutor*>&);
	ExchangeResult Begin(const ExchangeEpoch& epoch);
	void Wait(ExchangeResult& result);

private:
	const PartitionPlan& m_Plan;
	HaloTransport& m_Transport;
	std::map<uint32_t, DomainExecutor*> m_Executors;
	bool m_Active;
	bool m_HaveEpoch;
	ExchangeEpoch m_LastEpoch;
	std::map<uint32_t, HostBuffer> m_SendBuffers, m_ReceiveBuffers;
};

class FDTDCoordinator
{
public:
	FDTDCoordinator(const PartitionPlan&,
	    HaloTransport&,
	    CollectiveContext&,
	    const std::vector<DomainExecutor*>&,
	    uint64_t initial_step = 0);
	void Bootstrap();
	void Step();
	// Existing Engine extensions run between BeginPhase and FinishPhase. Halos
	// carry corrected fields only after every post/apply hook has completed.
	void BeginPhase(FieldPhase);
	void FinishPhase(FieldPhase);

	uint64_t NextStep() const
	{
		return m_Step;
	}

	bool Failed() const
	{
		return m_Failed;
	}

private:
	void Update(FieldPhase);
	const PartitionPlan& m_Plan;
	HaloTransport& m_Transport;
	CollectiveContext& m_Collectives;
	std::vector<DomainExecutor*> m_Executors;
	HaloExchange m_Exchange;
	uint64_t m_Step;
	bool m_Bootstrapped, m_Failed;
	bool m_PhasePending = false;
	FieldPhase m_NextPhase = FieldPhase::Voltage;
};

// Lifecycle always exists, even for a single domain. It is not a field relay.
class RunManager
{
public:
	RunManager(const PartitionPlan&, const std::vector<Placement>&, FDTDCoordinator& coordinator);
	void Start();
	void Advance(uint64_t steps);
	void Finish();
	void Run(uint64_t steps);

	bool Finished() const
	{
		return m_Finished;
	}

private:
	FDTDCoordinator& m_Coordinator;
	bool m_Started, m_Finished;
};

}
}
#endif

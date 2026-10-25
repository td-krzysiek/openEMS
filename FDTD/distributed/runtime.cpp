// SPDX-License-Identifier: GPL-3.0-or-later
#include "runtime.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <thread>

namespace openems
{
namespace distributed
{

namespace
{
class ReadyCompletion : public Completion
{
public:
	bool Ready() override
	{
		return true;
	}

	void Wait() override
	{
	}
};

struct LocalState
{
	bool complete = false;
	std::string error;
	HostBuffer send, receive;
};

class LocalCompletion : public Completion
{
public:
	LocalCompletion(std::shared_ptr<LocalState> state, std::function<void()> progress)
	    : m_State(state), m_Progress(progress)
	{
	}

	bool Ready() override
	{
		if (!m_State->error.empty())
		{
			throw std::runtime_error(m_State->error);
		}
		if (m_State->complete)
		{
			return true;
		}
		m_Progress();
		if (!m_State->error.empty())
		{
			throw std::runtime_error(m_State->error);
		}
		return m_State->complete;
	}

	void Wait() override
	{
		const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(30);
		while (!Ready())
		{
			if (std::chrono::steady_clock::now() > deadline)
			{
				throw std::runtime_error("local halo timeout: unmatched peer or failed progress");
			}
			std::this_thread::yield();
		}
	}

private:
	std::shared_ptr<LocalState> m_State;
	std::function<void()> m_Progress;
};

class UnpackCompletion : public Completion
{
public:
	UnpackCompletion(
	    CompletionPtr received, DomainExecutor* executor, HaloRecipe recipe, HostBuffer buffer)
	    : m_Received(received), m_Executor(executor), m_Recipe(recipe), m_Buffer(buffer)
	{
	}

	bool Ready() override
	{
		if (!m_Unpacked)
		{
			if (!m_Received->Ready())
			{
				return false;
			}
			m_Unpacked = m_Executor->Unpack(m_Recipe, m_Buffer);
			if (!m_Unpacked)
			{
				throw std::runtime_error("missing unpack completion");
			}
		}
		return m_Unpacked->Ready();
	}

	void Wait() override
	{
		m_Received->Wait();
		if (!m_Unpacked)
		{
			m_Unpacked = m_Executor->Unpack(m_Recipe, m_Buffer);
		}
		if (!m_Unpacked)
		{
			throw std::runtime_error("missing unpack completion");
		}
		m_Unpacked->Wait();
	}

private:
	CompletionPtr m_Received, m_Unpacked;
	DomainExecutor* m_Executor;
	HaloRecipe m_Recipe;
	HostBuffer m_Buffer;
};
}

CompletionPtr Completed()
{
	return std::make_shared<ReadyCompletion>();
}

struct LocalHaloTransport::Impl
{
	std::map<MessageId, std::shared_ptr<LocalState>> pending;
	unsigned int polls, delay;
	std::string failure;

	explicit Impl(unsigned int delay_) : polls(0), delay(delay_)
	{
	}

	void Progress()
	{
		if (!failure.empty())
		{
			throw std::runtime_error(failure);
		}
		if (++polls <= delay)
		{
			return;
		}
		for (auto it = pending.begin(); it != pending.end();)
		{
			auto state = it->second;
			if (!state->send || !state->receive)
			{
				++it;
				continue;
			}
			if (state->send->size() != state->receive->size())
			{
				state->error = "local halo sample count mismatch";
				throw std::runtime_error(state->error);
			}
			std::copy(state->send->begin(), state->send->end(), state->receive->begin());
			state->complete = true;
			state->send.reset();
			state->receive.reset();
			it = pending.erase(it);
		}
	}

	void Fail(const std::string& reason)
	{
		failure = reason.empty() ? "local halo transport failed" : reason;
		for (auto& item : pending)
		{
			item.second->error = failure;
		}
	}
};

LocalHaloTransport::LocalHaloTransport(unsigned int delay) : m_Impl(std::make_shared<Impl>(delay))
{
}

LocalHaloTransport::~LocalHaloTransport()
{
	m_Impl->Fail("halo transport destroyed with outstanding work");
}

CompletionPtr LocalHaloTransport::Receive(const MessageId& id, HostBuffer buffer)
{
	if (!buffer || !m_Impl->failure.empty())
	{
		throw std::runtime_error("invalid receive or failed local transport");
	}
	auto& state = m_Impl->pending[id];
	if (!state)
	{
		state = std::make_shared<LocalState>();
	}
	if (state->receive)
	{
		throw std::logic_error("duplicate local receive");
	}
	state->receive = buffer;
	std::weak_ptr<Impl> impl = m_Impl;
	return std::make_shared<LocalCompletion>(state,
	    [impl]()
	{
		auto p = impl.lock();
		if (!p)
		{
			throw std::runtime_error("local transport no longer exists");
		}
		p->Progress();
	});
}

CompletionPtr LocalHaloTransport::Send(const MessageId& id, HostBuffer buffer)
{
	if (!buffer || !m_Impl->failure.empty())
	{
		throw std::runtime_error("invalid send or failed local transport");
	}
	auto& state = m_Impl->pending[id];
	if (!state)
	{
		state = std::make_shared<LocalState>();
	}
	if (state->send)
	{
		throw std::logic_error("duplicate local send");
	}
	state->send = buffer;
	std::weak_ptr<Impl> impl = m_Impl;
	return std::make_shared<LocalCompletion>(state,
	    [impl]()
	{
		auto p = impl.lock();
		if (!p)
		{
			throw std::runtime_error("local transport no longer exists");
		}
		p->Progress();
	});
}

void LocalHaloTransport::Progress()
{
	m_Impl->Progress();
}

void LocalHaloTransport::Fail(const std::string& reason)
{
	m_Impl->Fail(reason);
}

double LocalCollectiveContext::Minimum(double value)
{
	if (!std::isfinite(value))
	{
		throw std::invalid_argument("nonfinite collective input");
	}
	return value;
}

double LocalCollectiveContext::Sum(double value)
{
	return Minimum(value);
}

void LocalCollectiveContext::Fail(const std::string& reason)
{
	throw std::runtime_error(reason);
}

HaloExchange::HaloExchange(const PartitionPlan& plan,
    HaloTransport& transport,
    const std::vector<DomainExecutor*>& executors)
    : m_Plan(plan), m_Transport(transport), m_Active(false), m_HaveEpoch(false),
      m_LastEpoch{ 0, 0, FieldPhase::Current }
{
	for (auto* e : executors)
	{
		if (!e || e->Id() >= plan.Domains().size() || !m_Executors.emplace(e->Id(), e).second)
		{
			throw std::invalid_argument("invalid or duplicate local executor");
		}
	}
	for (const auto& recipe : plan.Halos())
	{
		if (m_Executors.count(recipe.source))
		{
			m_SendBuffers[recipe.channel] = std::make_shared<std::vector<float>>(recipe.Samples());
		}
		if (m_Executors.count(recipe.destination))
		{
			m_ReceiveBuffers[recipe.channel] =
			    std::make_shared<std::vector<float>>(recipe.Samples());
		}
	}
}

ExchangeResult HaloExchange::Begin(const ExchangeEpoch& epoch)
{
	if (m_Active || epoch.plan != m_Plan.Id())
	{
		throw std::logic_error("overlapping exchange or wrong plan ID");
	}
	if (m_HaveEpoch)
	{
		const bool after_v = m_LastEpoch.phase == FieldPhase::Voltage &&
		    epoch.phase == FieldPhase::Current && epoch.step == m_LastEpoch.step;
		const bool after_i = m_LastEpoch.phase == FieldPhase::Current &&
		    epoch.phase == FieldPhase::Voltage && m_LastEpoch.step != UINT64_MAX &&
		    epoch.step == m_LastEpoch.step + 1;
		if (!after_v && !after_i)
		{
			throw std::logic_error("stale, skipped or wrong-phase exchange epoch");
		}
	}
	else if (epoch.phase != FieldPhase::Current)
	{
		throw std::logic_error("bootstrap current halos before voltage update");
	}
	m_Active = true;
	m_HaveEpoch = true;
	m_LastEpoch = epoch;
	ExchangeResult result;
	result.m_Epoch = epoch;
	try
	{
		// All receives are posted before any send wait, including multi-rank chains.
		for (const auto& recipe : m_Plan.Halos())
		{
			if (recipe.phase == epoch.phase && m_Executors.count(recipe.destination))
			{
				MessageId id = { epoch, recipe.channel, recipe.source, recipe.destination };
				auto b = m_ReceiveBuffers.at(recipe.channel);
				auto received = m_Transport.Receive(id, b);
				if (!received)
				{
					throw std::runtime_error("missing receive completion");
				}
				result.m_ReceiveReady.push_back(std::make_shared<UnpackCompletion>(
				    received, m_Executors.at(recipe.destination), recipe, b));
			}
		}
		for (const auto& recipe : m_Plan.Halos())
		{
			if (recipe.phase == epoch.phase && m_Executors.count(recipe.source))
			{
				MessageId id = { epoch, recipe.channel, recipe.source, recipe.destination };
				auto b = m_SendBuffers.at(recipe.channel);
				auto packed = m_Executors.at(recipe.source)->Pack(recipe, b);
				if (!packed)
				{
					throw std::runtime_error("missing pack completion");
				}
				packed->Wait();
				auto sent = m_Transport.Send(id, b);
				if (!sent)
				{
					throw std::runtime_error("missing send completion");
				}
				result.m_SendRelease.push_back(sent);
			}
		}
	}
	catch (const std::exception& error)
	{
		m_Transport.Fail(error.what());
		throw;
	}
	return result;
}

void HaloExchange::Wait(ExchangeResult& result)
{
	if (!m_Active || !(result.m_Epoch == m_LastEpoch))
	{
		throw std::logic_error("no matching active halo exchange");
	}
	try
	{
		for (auto& token : result.m_ReceiveReady)
		{
			token->Wait();
		}
		for (auto& token : result.m_SendRelease)
		{
			token->Wait();
		}
		m_Transport.Progress();
		m_Active = false;
	}
	catch (const std::exception& error)
	{
		m_Transport.Fail(error.what());
		throw;
	}
}

FDTDCoordinator::FDTDCoordinator(const PartitionPlan& plan,
    HaloTransport& transport,
    CollectiveContext& collectives,
    const std::vector<DomainExecutor*>& executors,
    uint64_t initial_step)
    : m_Plan(plan), m_Transport(transport), m_Collectives(collectives), m_Executors(executors),
      m_Exchange(plan, transport, executors), m_Step(initial_step), m_Bootstrapped(false),
      m_Failed(false)
{
}

void FDTDCoordinator::Bootstrap()
{
	if (m_Failed || m_Bootstrapped)
	{
		throw std::logic_error("invalid coordinator bootstrap");
	}
	try
	{
		for (const auto& domain : m_Plan.Domains())
		{
			const size_t local = std::count_if(m_Executors.begin(),
			    m_Executors.end(),
			    [&](DomainExecutor* e)
			{
				return e->Id() == domain.id;
			});
			if (m_Collectives.Sum(double(local)) != 1.)
			{
				throw std::logic_error(
				    "each domain must have exactly one active executor globally");
			}
		}
		double local_dt = std::numeric_limits<double>::max();
		for (auto* e : m_Executors)
		{
			local_dt = std::min(local_dt, e->Timestep());
		}
		const double global_dt = m_Collectives.Minimum(local_dt);
		bool same = std::isfinite(global_dt) && global_dt > 0;
		for (auto* e : m_Executors)
		{
			same = same && e->Timestep() == global_dt;
		}
		if (!m_Collectives.All(same))
		{
			throw std::logic_error("executors disagree on the global operator timestep");
		}
		auto exchange = m_Exchange.Begin({ m_Plan.Id(), m_Step, FieldPhase::Current });
		m_Exchange.Wait(exchange);
		m_Bootstrapped = true;
	}
	catch (const std::exception& error)
	{
		m_Failed = true;
		m_Transport.Fail(error.what());
		m_Collectives.Fail(error.what());
		throw;
	}
}

void FDTDCoordinator::Update(FieldPhase phase)
{
	std::vector<CompletionPtr> work;
	for (auto* executor : m_Executors)
	{
		auto token = executor->Update(phase, m_Step);
		if (!token)
		{
			throw std::runtime_error("missing execution completion");
		}
		work.push_back(token);
	}
	for (auto& token : work)
	{
		token->Wait();
	}
}

void FDTDCoordinator::Step()
{
	BeginPhase(FieldPhase::Voltage);
	FinishPhase(FieldPhase::Voltage);
	BeginPhase(FieldPhase::Current);
	FinishPhase(FieldPhase::Current);
}

void FDTDCoordinator::BeginPhase(FieldPhase phase)
{
	if (!m_Bootstrapped || m_Failed || m_Step == UINT64_MAX)
	{
		throw std::logic_error("coordinator not bootstrapped, failed or timestep overflow");
	}
	if (m_PhasePending || phase != m_NextPhase)
	{
		throw std::logic_error("overlapping or out-of-order domain half-step");
	}
	try
	{
		Update(phase);
		m_PhasePending = true;
	}
	catch (const std::exception& error)
	{
		m_Failed = true;
		m_Transport.Fail(error.what());
		m_Collectives.Fail(error.what());
		throw;
	}
}

void FDTDCoordinator::FinishPhase(FieldPhase phase)
{
	if (!m_PhasePending || m_Failed || phase != m_NextPhase)
	{
		throw std::logic_error("no matching domain half-step to finish");
	}
	try
	{
		auto exchange = m_Exchange.Begin({ m_Plan.Id(), m_Step + 1, phase });
		m_Exchange.Wait(exchange);
		for (auto* executor : m_Executors)
		{
			executor->Observe(phase, m_Step);
		}
		m_PhasePending = false;
		m_NextPhase = phase == FieldPhase::Voltage ? FieldPhase::Current : FieldPhase::Voltage;
		if (phase == FieldPhase::Current)
		{
			++m_Step;
		}
	}
	catch (const std::exception& error)
	{
		m_Failed = true;
		m_Transport.Fail(error.what());
		m_Collectives.Fail(error.what());
		throw;
	}
}

RunManager::RunManager(const PartitionPlan& plan,
    const std::vector<Placement>& placement,
    FDTDCoordinator& coordinator)
    : m_Coordinator(coordinator), m_Started(false), m_Finished(false)
{
	ValidatePlacement(plan, placement);
}

void RunManager::Start()
{
	if (m_Started)
	{
		throw std::logic_error("run already started");
	}
	m_Coordinator.Bootstrap();
	m_Started = true;
}

void RunManager::Advance(uint64_t steps)
{
	if (!m_Started || m_Finished || !steps || steps > UINT64_MAX - m_Coordinator.NextStep())
	{
		throw std::logic_error("run inactive, finished or invalid step count");
	}
	for (uint64_t n = 0; n < steps; ++n)
	{
		m_Coordinator.Step();
	}
}

void RunManager::Finish()
{
	if (!m_Started || m_Finished || m_Coordinator.Failed())
	{
		throw std::logic_error("run inactive, failed or already finished");
	}
	m_Finished = true;
}

void RunManager::Run(uint64_t steps)
{
	if (!steps || steps > UINT64_MAX - m_Coordinator.NextStep())
	{
		throw std::invalid_argument("invalid step count");
	}
	Start();
	Advance(steps);
	Finish();
}

}
}

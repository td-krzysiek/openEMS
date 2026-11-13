// SPDX-License-Identifier: GPL-3.0-or-later
#ifdef OPENEMS_WITH_DOMAIN_RUNTIME
#include "gpu_backend_domains.h"
#include "distributed/runtime.h"
#include "domain_performance.h"
#include "extensions/engine_ext_excitation.h"
#include "extensions/engine_ext_lorentzmaterial.h"
#include "extensions/engine_ext_lumpedRLC.h"
#include "extensions/engine_ext_opaquesurface.h"
#include "extensions/engine_ext_upml.h"
#include "gpu_backend_reference.h"
#include "operator.h"
#include "operator_cylinder.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

using namespace openems::distributed;

namespace
{
class BackendCompletion : public Completion
{
public:
	BackendCompletion(GPU_Backend* backend, DomainPerformance* profile, size_t domain,
		DomainPerformance::HostMetric metric = DomainPerformance::UpdateWait)
		: m_Backend(backend), m_Profile(profile), m_Domain(domain), m_Metric(metric), m_Done(false)
	{
	}
	bool Ready() override
	{
		return m_Done;
	}
	void Wait() override
	{
		if (!m_Done)
		{
			m_Backend->Activate();
			const bool sampled = m_Profile && m_Profile->Active();
			const auto started = sampled ? DomainPerformance::Now() : DomainPerformance::Time();
			m_Backend->Synchronize();
			if (sampled)
			{
				m_Profile->Host(m_Domain, m_Metric, started);
			}
			m_Done = true;
		}
	}

private:
	GPU_Backend* m_Backend;
	DomainPerformance* m_Profile;
	size_t m_Domain;
	DomainPerformance::HostMetric m_Metric;
	bool m_Done;
};

// Buffer ownership remains with the executor until the last operation drains.
class HaloCompletion : public Completion
{
public:
	HaloCompletion(GPU_HaloBuffer* native, HostBuffer host, bool download)
	    : m_Native(native), m_Host(host), m_Download(download)
	{
	}
	bool Ready() override
	{
		if (!m_Native->Ready())
		{
			return false;
		}
		Finish();
		return true;
	}
	void Wait() override
	{
		m_Native->Wait();
		Finish();
	}
private:
	void Finish()
	{
		if (!m_Done && m_Download && m_Host->data() != m_Native->Data())
		{
			std::copy(m_Native->Data(), m_Native->Data() + m_Host->size(), m_Host->begin());
		}
		m_Done = true;
	}
	GPU_HaloBuffer* m_Native;
	HostBuffer m_Host;
	bool m_Download, m_Done = false;
};

class BackendExecutor : public DomainExecutor
{
public:
	BackendExecutor(Domain domain, GPU_Backend* backend, double timestep, DomainPerformance* profile,
	    const std::vector<HaloRecipe>& recipes, bool async, bool borrowed)
		: m_Domain(domain), m_Backend(backend), m_Profile(profile), m_Timestep(timestep),
		  m_Borrowed(borrowed)
	{
		if (async)
		{
			for (const auto& recipe : recipes)
			{
				if (recipe.source == Id() || recipe.destination == Id())
				{
					std::shared_ptr<GPU_HaloBuffer> buffer(backend->NewHaloBuffer(recipe.Samples()));
					if (!buffer || (borrowed && !buffer->CanUploadBorrowed()))
					{
						throw std::runtime_error("selected backend has no asynchronous halo adapter");
					}
					m_Buffers.emplace(recipe.channel, std::move(buffer));
				}
			}
		}
	}
	bool CanOverlapCurrentInterior() const override
	{
		return !m_Buffers.empty() && m_Backend->CanOverlapCurrentInterior();
	}
	CompletionPtr UpdateCurrentRegion(bool interior, uint64_t) override
	{
		m_Backend->Activate();
		auto* timer = m_Profile ? m_Profile->Begin(Id(), false) : nullptr;
		m_Backend->UpdateCurrentRegion(interior);
		if (m_Profile)
		{
			m_Profile->End(Id(), timer);
		}
		return std::make_shared<BackendCompletion>(m_Backend, m_Profile, Id());
	}
	uint32_t Id() const override
	{
		return m_Domain.id;
	}
	double Timestep() const override
	{
		return m_Timestep;
	}
	CompletionPtr Update(FieldPhase phase, uint64_t) override
	{
		m_Backend->Activate();
		auto* timing = m_Profile ? m_Profile->Begin(Id(), false) : nullptr;
		if (phase == FieldPhase::Voltage)
		{
			m_Backend->UpdateVoltages();
		}
		else
		{
			m_Backend->UpdateCurrents();
		}
		if (m_Profile)
		{
			m_Profile->End(Id(), timing);
		}
		return std::make_shared<BackendCompletion>(m_Backend, m_Profile, Id());
	}
	HostBuffer CreateHaloBuffer(const HaloRecipe &recipe, bool sending) override
	{
		if (m_Borrowed && sending)
		{
			auto buffer = m_Buffers.at(recipe.channel);
			return std::make_shared<HostSamples>(recipe.Samples(), buffer->Data(), buffer);
		}
		auto values = std::make_shared<HostSamples>(recipe.Samples());
		if (m_Borrowed)
		{
			values->AllowBorrowing();
		}
		return values;
	}
	CompletionPtr Pack(const HaloRecipe& recipe, HostBuffer values) override
	{
		const bool sampled = m_Profile && m_Profile->Active();
		const auto started = sampled ? DomainPerformance::Now() : DomainPerformance::Time();
		m_Backend->Activate();
		if (!m_Buffers.empty())
		{
			auto* buffer = m_Buffers.at(recipe.channel).get();
			const size_t offsets[2] = { m_Domain.LocalIndex(1, recipe.source_x, 0, 0),
				m_Domain.LocalIndex(2, recipe.source_x, 0, 0) };
			buffer->Download(recipe.phase == FieldPhase::Current, offsets, recipe.plane_samples);
			if (sampled)
			{
				m_Profile->Host(Id(), DomainPerformance::Pack, started, values->size() * sizeof(float));
			}
			return std::make_shared<HaloCompletion>(buffer, values, true);
		}
		for (unsigned int c = 1; c < 3; ++c)
		{
			const size_t offset = m_Domain.LocalIndex(c, recipe.source_x, 0, 0);
			if (!m_Backend->DownloadRange(recipe.phase == FieldPhase::Current,
			        offset,
			        recipe.plane_samples,
			        values->data() + (c - 1) * recipe.plane_samples))
			{
				throw std::runtime_error("domain executor cannot pack a halo");
			}
		}
		if (sampled)
		{
			m_Profile->Host(Id(), DomainPerformance::Pack, started, values->size() * sizeof(float));
		}
		return Completed();
	}
	CompletionPtr Unpack(const HaloRecipe& recipe, HostBuffer values) override
	{
		const bool sampled = m_Profile && m_Profile->Active();
		const auto started = sampled ? DomainPerformance::Now() : DomainPerformance::Time();
		m_Backend->Activate();
		if (!m_Buffers.empty())
		{
			auto* buffer = m_Buffers.at(recipe.channel).get();
			const size_t offsets[2] = { m_Domain.LocalIndex(1, recipe.destination_x, 0, 0),
				m_Domain.LocalIndex(2, recipe.destination_x, 0, 0) };
			if (values->Owner())
			{
				buffer->UploadBorrowed(recipe.phase == FieldPhase::Current, offsets, recipe.plane_samples,
									   values->data(), values->Owner());
			}
			else
			{
				std::copy(values->begin(), values->end(), buffer->Data());
				buffer->Upload(recipe.phase == FieldPhase::Current, offsets, recipe.plane_samples);
			}
			if (sampled)
			{
				m_Profile->Host(Id(), DomainPerformance::Unpack, started, values->size() * sizeof(float));
			}
			return std::make_shared<HaloCompletion>(buffer, values, false);
		}
		for (unsigned int c = 1; c < 3; ++c)
		{
			const size_t offset = m_Domain.LocalIndex(c, recipe.destination_x, 0, 0);
			if (!m_Backend->UploadRange(recipe.phase == FieldPhase::Current,
			        offset,
			        recipe.plane_samples,
			        values->data() + (c - 1) * recipe.plane_samples))
			{
				throw std::runtime_error("domain executor cannot unpack a halo");
			}
		}
		if (sampled)
		{
			m_Profile->Host(
				Id(), DomainPerformance::Unpack, started, values->size() * sizeof(float));
		}
		return std::make_shared<BackendCompletion>(
			m_Backend, m_Profile, Id(), DomainPerformance::UnpackWait);
	}

private:
	Domain m_Domain;
	GPU_Backend* m_Backend;
	DomainPerformance* m_Profile;
	double m_Timestep;
	bool m_Borrowed;
	std::map<uint32_t, std::shared_ptr<GPU_HaloBuffer>> m_Buffers;
};

class DomainExtension : public GPU_Extension
{
public:
	DomainExtension(DomainPerformance* profile, const std::string& name)
		: m_Profile(profile), m_ProfileId(profile ? profile->RegisterExtension(name) : DomainPerformance::NoExtension)
	{
	}
	DomainPerformance* m_Profile;
	size_t m_ProfileId;
	std::vector<std::pair<GPU_Backend*, std::unique_ptr<GPU_Extension>>> children;
	~DomainExtension() override
	{
		for (auto& child : children)
		{
			try
			{
				child.first->Activate();
			}
			catch (...)
			{
			}
			child.second.reset();
		}
	}
	void Invoke(void (GPU_Extension::*hook)(), DomainPerformance::ExtensionHook label)
	{
		for (size_t i = 0; i < children.size(); ++i)
		{
			auto& child = children[i];
			child.first->Activate();
			auto* timing = m_Profile ? m_Profile->BeginExtension(i, m_ProfileId, label) : nullptr;
			(child.second.get()->*hook)();
			if (m_Profile)
			{
				m_Profile->End(i, timing);
			}
		}
	}
	void DoPreVoltageUpdates() override
	{
		Invoke(&GPU_Extension::DoPreVoltageUpdates, DomainPerformance::PreVoltage);
	}
	void DoPostVoltageUpdates() override
	{
		Invoke(&GPU_Extension::DoPostVoltageUpdates, DomainPerformance::PostVoltage);
	}
	void Apply2Voltages() override
	{
		Invoke(&GPU_Extension::Apply2Voltages, DomainPerformance::ApplyVoltage);
	}
	void DoPreCurrentUpdates() override
	{
		Invoke(&GPU_Extension::DoPreCurrentUpdates, DomainPerformance::PreCurrent);
	}
	void DoPostCurrentUpdates() override
	{
		Invoke(&GPU_Extension::DoPostCurrentUpdates, DomainPerformance::PostCurrent);
	}
	void Apply2Current() override
	{
		Invoke(&GPU_Extension::Apply2Current, DomainPerformance::ApplyCurrent);
	}
	void Synchronize() override
	{
		Invoke(&GPU_Extension::Synchronize, DomainPerformance::SynchronizeHook);
	}
};
}

struct GPU_Backend_Domains::Impl
{
	std::vector<std::unique_ptr<GPU_Backend>> devices;
	std::vector<uint32_t> cuts;
	std::unique_ptr<PartitionPlan> plan;
	std::unique_ptr<DomainPerformance> profile;
	LocalHaloTransport transport;
	LocalCollectiveContext collectives;
	std::vector<std::unique_ptr<BackendExecutor>> executors;
	std::unique_ptr<FDTDCoordinator> coordinator;
	std::vector<float> lower_boundary;
	void MirrorLowerCurrent()
	{
		auto* backend = devices.front().get();
		const auto& domain = plan->Domains().front();
		backend->Activate();
		for (unsigned int c = 1; c < 3; ++c)
		{
			const size_t source = domain.LocalIndex(c, 1, 0, 0);
			const size_t target = domain.LocalIndex(c, 0, 0, 0);
			const size_t samples = domain.grid.Plane();
			if (backend->CopyRange(true, source, target, samples))
			{
				continue;
			}
			lower_boundary.resize(samples);
			if (!backend->DownloadRange(true, source, samples, lower_boundary.data()) ||
			    !backend->UploadRange(true, target, samples, lower_boundary.data()))
			{
				throw std::runtime_error("cannot apply physical lower boundary");
			}
		}
		backend->Synchronize();
	}
};

GPU_Backend_Domains::GPU_Backend_Domains(
    std::vector<std::unique_ptr<GPU_Backend>> devices, std::vector<uint32_t> cuts)
    : d(new Impl())
{
	if (devices.size() < 2)
	{
		throw std::invalid_argument("select at least two domain executors");
	}
	if (!cuts.empty() && cuts.size() != devices.size() + 1)
	{
		throw std::invalid_argument("domain cuts must contain one more boundary than devices");
	}
	d->cuts = std::move(cuts);
	d->devices = std::move(devices);
}

GPU_Backend_Domains::~GPU_Backend_Domains()
{
	// Timers belong to child contexts and must be destroyed before those contexts.
	d->coordinator.reset();
	d->transport.Fail("domain backend destroyed");
	d->executors.clear();
	d->profile.reset();
	for (auto& backend : d->devices)
	{
		try
		{
			backend->Activate();
		}
		catch (...)
		{
		}
		backend.reset();
	}
}

std::string GPU_Backend_Domains::GetName() const
{
	std::string name = "domains [";
	for (const auto& backend : d->devices)
	{
		name += backend->GetName() + "; ";
	}
	return name + "]";
}

bool GPU_Backend_Domains::Init(const Operator* op)
{
	if (dynamic_cast<const Operator_Cylinder*>(op))
	{
		throw std::invalid_argument("domain backend requires a Cartesian operator");
	}
	GridShape shape = {
		op->GetNumberOfLines(0, true), op->GetNumberOfLines(1, true), op->GetNumberOfLines(2, true)
	};
	std::vector<uint32_t> cuts = d->cuts;
	if (cuts.empty())
	{
		for (size_t i = 0; i <= d->devices.size(); ++i)
		{
			cuts.push_back(uint32_t(uint64_t(shape.nx) * i / d->devices.size()));
		}
	}
	const char* setting = std::getenv("OPENEMS_DOMAIN_HALOS");
	const std::string mode = setting ? setting : "serial";
	if (mode != "serial" && mode != "async" && mode != "overlap")
	{
		throw std::invalid_argument("OPENEMS_DOMAIN_HALOS must be serial, async or overlap");
	}
	const bool async = mode != "serial";
	const char *host_setting = std::getenv("OPENEMS_DOMAIN_HOST_BUFFERS");
	const std::string host_mode = host_setting ? host_setting : "copy";
	if ((host_mode != "copy" && host_mode != "borrowed") || (host_mode == "borrowed" && !async))
	{
		throw std::invalid_argument("borrowed domain host buffers require async or overlap halos");
	}
	const bool borrowed = host_mode == "borrowed";
	std::cout << "Domain host buffers: " << host_mode << std::endl;
	std::cout << "Domain halo mode: " << mode << std::endl;
	d->plan.reset(new PartitionPlan(1, shape, cuts));
	std::vector<GPU_Backend*> profile_devices;
	for (const auto& device : d->devices)
	{
		profile_devices.push_back(device.get());
	}
	d->profile = DomainPerformance::FromEnvironment(*d->plan, profile_devices);
	std::vector<DomainExecutor*> executors;
	for (const auto& domain : d->plan->Domains())
	{
		auto* backend = d->devices[domain.id].get();
		backend->Activate();
		if (!backend->InitDomain(op, domain.begin, domain.end))
		{
			throw std::runtime_error("selected executor does not support domain initialization");
		}
		d->executors.emplace_back(new BackendExecutor(domain, backend, op->GetTimestep(), d->profile.get(),
													  d->plan->Halos(), async, borrowed));
		executors.push_back(d->executors.back().get());
		std::cout << "Domain " << domain.id << ": owns X[" << domain.begin << "," << domain.end
		          << "), storage nodes " << domain.StorageNodes() << ", " << backend->GetName()
		          << std::endl;
	}
	d->coordinator.reset(new FDTDCoordinator(*d->plan, d->transport, d->collectives, executors, 0,
	        mode == "overlap"));
	d->coordinator->Bootstrap();
	d->MirrorLowerCurrent();
	if (d->profile)
	{
		d->profile->StartStep(d->coordinator->NextStep());
	}
	return true;
}

void GPU_Backend_Domains::UpdateVoltages()
{
	d->coordinator->BeginPhase(FieldPhase::Voltage);
}
void GPU_Backend_Domains::UpdateCurrents()
{
	d->coordinator->BeginPhase(FieldPhase::Current);
	if (d->profile)
	{
		d->profile->ExchangeSeconds(d->coordinator->LastDeferredExchangeSeconds());
	}
}
void GPU_Backend_Domains::FinishVoltagePhase()
{
	const bool sampled = d->profile && d->profile->Active();
	const auto started = sampled ? DomainPerformance::Now() : DomainPerformance::Time();
	d->coordinator->FinishPhase(FieldPhase::Voltage);
	if (sampled)
	{
		d->profile->Exchange(started);
	}
}
void GPU_Backend_Domains::FinishCurrentPhase()
{
	const bool sampled = d->profile && d->profile->Active();
	const auto started = sampled ? DomainPerformance::Now() : DomainPerformance::Time();
	d->coordinator->FinishPhase(FieldPhase::Current);
	if (sampled)
	{
		d->profile->Exchange(started);
	}
	const auto boundary = sampled ? DomainPerformance::Now() : DomainPerformance::Time();
	d->MirrorLowerCurrent();
	if (sampled)
	{
		d->profile->Host(0, DomainPerformance::Boundary, boundary);
	}
	if (d->profile)
	{
		d->profile->FinishStep();
		d->profile->StartStep(d->coordinator->NextStep());
	}
}

bool GPU_Backend_Domains::DownloadRange(bool currents, size_t offset, size_t count, float* dst)
{
	const auto& shape = d->plan->Shape();
	const size_t nodes = shape.Nodes(), plane = shape.Plane();
	if (offset > 3 * nodes || count > 3 * nodes - offset)
	{
		throw std::out_of_range("global domain download range");
	}
	while (count)
	{
		const unsigned int c = unsigned(offset / nodes);
		const size_t node = offset % nodes;
		const unsigned int x = unsigned(node / plane);
		const auto& domain = d->plan->Domains()[d->plan->Owner(x)];
		const size_t n = std::min(count, size_t(domain.end) * plane - node);
		const size_t local = domain.LocalIndex(c, domain.ToLocal(x), 0, 0) + node % plane;
		auto* backend = d->devices[domain.id].get();
		backend->Activate();
		if (!backend->DownloadRange(currents, local, n, dst))
		{
			return false;
		}
		offset += n;
		count -= n;
		dst += n;
	}
	return true;
}

void GPU_Backend_Domains::Download(bool currents, ArrayLib::ArrayNIJK<float>& field)
{
	if (!DownloadRange(currents, 0, field.size(), field.data()))
	{
		throw std::runtime_error("cannot gather domain fields");
	}
}

void GPU_Backend_Domains::Upload(bool currents, const ArrayLib::ArrayNIJK<float>& field)
{
	const size_t plane = d->plan->Shape().Plane(), nodes = d->plan->Shape().Nodes();
	for (const auto& domain : d->plan->Domains())
	{
		auto* backend = d->devices[domain.id].get();
		backend->Activate();
		for (unsigned int c = 0; c < 3; ++c)
		{
			for (unsigned int x = 0; x < domain.end - domain.begin + 2; ++x)
			{
				const int gx =
				    std::max(0, std::min(int(domain.grid.nx) - 1, int(domain.begin) + int(x) - 1));
				if (!backend->UploadRange(currents,
				        domain.LocalIndex(c, x, 0, 0),
				        plane,
				        field.data() + c * nodes + size_t(gx) * plane))
				{
					throw std::runtime_error("cannot scatter domain fields");
				}
			}
		}
	}
	Synchronize();
}

void GPU_Backend_Domains::DownloadVoltages(ArrayLib::ArrayNIJK<float>& field)
{
	Download(false, field);
}
void GPU_Backend_Domains::DownloadCurrents(ArrayLib::ArrayNIJK<float>& field)
{
	Download(true, field);
}
void GPU_Backend_Domains::UploadVoltages(const ArrayLib::ArrayNIJK<float>& field)
{
	Upload(false, field);
}
void GPU_Backend_Domains::UploadCurrents(const ArrayLib::ArrayNIJK<float>& field)
{
	Upload(true, field);
}

void GPU_Backend_Domains::Synchronize()
{
	for (auto& backend : d->devices)
	{
		backend->Activate();
		backend->Synchronize();
	}
}

GPU_Extension* GPU_Backend_Domains::CreateExtension(Engine_Extension* ext, Engine* engine)
{
	std::unique_ptr<DomainExtension> result(new DomainExtension(d->profile.get(),
		d->profile ? ext->GetExtensionName() : std::string()));
	for (auto& backend : d->devices)
	{
		backend->Activate();
		if (!dynamic_cast<GPU_Backend_Reference*>(backend.get()) &&
		    !dynamic_cast<Engine_Ext_UPML*>(ext) && !dynamic_cast<Engine_Ext_Excitation*>(ext) &&
		    !dynamic_cast<Engine_Ext_LumpedRLC*>(ext) &&
		    !dynamic_cast<Engine_Ext_OpaqueSurface*>(ext) &&
		    !dynamic_cast<Engine_Ext_LorentzMaterial*>(ext))
		{
			throw std::runtime_error("unsupported extension for native device-domain execution");
		}
		std::unique_ptr<GPU_Extension> child(backend->CreateExtension(ext, engine));
		if (!child)
		{
			if (!dynamic_cast<GPU_Backend_Reference*>(backend.get()))
			{
				throw std::runtime_error(
				    "extension has no qualified device-domain implementation; host fallback is disabled for multi-GPU execution");
			}
			return nullptr;
		}
		result->children.emplace_back(backend.get(), std::move(child));
	}
	return result.release();
}

GPU_Backend* GPU_Backend_Domains::NewSubGridBackend()
{
	throw std::runtime_error("cylindrical multi-grid is unsupported by the domain backend");
}

bool GPU_Backend_Domains::CalcFastEnergy(
    const unsigned int numNodes[3], double& electric, double& magnetic)
{
	electric = magnetic = 0;
	for (auto& backend : d->devices)
	{
		backend->Activate();
		double e = 0, h = 0;
		if (!backend->CalcFastEnergy(numNodes, e, h))
		{
			return false;
		}
		electric += e;
		magnetic += h;
	}
	return true;
}
#endif

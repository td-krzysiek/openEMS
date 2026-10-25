// SPDX-License-Identifier: GPL-3.0-or-later
#ifdef OPENEMS_WITH_DOMAIN_RUNTIME
#include "gpu_backend_domains.h"
#include "distributed/runtime.h"
#include "operator.h"
#include "operator_cylinder.h"
#include "gpu_backend_reference.h"
#include <iostream>
#include <stdexcept>

using namespace openems::distributed;

namespace
{
class BackendCompletion : public Completion
{
public:
	explicit BackendCompletion(GPU_Backend* backend) : m_Backend(backend), m_Done(false)
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
			m_Backend->Synchronize();
			m_Done = true;
		}
	}

private:
	GPU_Backend* m_Backend;
	bool m_Done;
};

class BackendExecutor : public DomainExecutor
{
public:
	BackendExecutor(Domain domain, GPU_Backend* backend, double timestep)
	    : m_Domain(domain), m_Backend(backend), m_Timestep(timestep)
	{
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
		if (phase == FieldPhase::Voltage)
		{
			m_Backend->UpdateVoltages();
		}
		else
		{
			m_Backend->UpdateCurrents();
		}
		return std::make_shared<BackendCompletion>(m_Backend);
	}
	CompletionPtr Pack(const HaloRecipe& recipe, HostBuffer values) override
	{
		m_Backend->Activate();
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
		return Completed();
	}
	CompletionPtr Unpack(const HaloRecipe& recipe, HostBuffer values) override
	{
		m_Backend->Activate();
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
		return std::make_shared<BackendCompletion>(m_Backend);
	}

private:
	Domain m_Domain;
	GPU_Backend* m_Backend;
	double m_Timestep;
};

class DomainExtension : public GPU_Extension
{
public:
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
	void Invoke(void (GPU_Extension::*hook)())
	{
		for (auto& child : children)
		{
			child.first->Activate();
			(child.second.get()->*hook)();
		}
	}
	void DoPreVoltageUpdates() override
	{
		Invoke(&GPU_Extension::DoPreVoltageUpdates);
	}
	void DoPostVoltageUpdates() override
	{
		Invoke(&GPU_Extension::DoPostVoltageUpdates);
	}
	void Apply2Voltages() override
	{
		Invoke(&GPU_Extension::Apply2Voltages);
	}
	void DoPreCurrentUpdates() override
	{
		Invoke(&GPU_Extension::DoPreCurrentUpdates);
	}
	void DoPostCurrentUpdates() override
	{
		Invoke(&GPU_Extension::DoPostCurrentUpdates);
	}
	void Apply2Current() override
	{
		Invoke(&GPU_Extension::Apply2Current);
	}
	void Synchronize() override
	{
		Invoke(&GPU_Extension::Synchronize);
	}
};
}

struct GPU_Backend_Domains::Impl
{
	std::vector<std::unique_ptr<GPU_Backend>> devices;
	std::unique_ptr<PartitionPlan> plan;
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
			if (!backend->DownloadRange(
			        true, source, lower_boundary.size(), lower_boundary.data()) ||
			    !backend->UploadRange(true, target, lower_boundary.size(), lower_boundary.data()))
			{
				throw std::runtime_error("cannot apply physical lower boundary");
			}
		}
		backend->Synchronize();
	}
};

GPU_Backend_Domains::GPU_Backend_Domains(std::vector<std::unique_ptr<GPU_Backend>> devices)
    : d(new Impl())
{
	if (devices.size() < 2)
	{
		throw std::invalid_argument("select at least two domain executors");
	}
	d->devices = std::move(devices);
}

GPU_Backend_Domains::~GPU_Backend_Domains()
{
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
	std::vector<uint32_t> cuts;
	for (size_t i = 0; i <= d->devices.size(); ++i)
	{
		cuts.push_back(uint32_t(uint64_t(shape.nx) * i / d->devices.size()));
	}
	d->plan.reset(new PartitionPlan(1, shape, cuts));
	std::vector<DomainExecutor*> executors;
	for (const auto& domain : d->plan->Domains())
	{
		auto* backend = d->devices[domain.id].get();
		backend->Activate();
		if (!backend->InitDomain(op, domain.begin, domain.end))
		{
			throw std::runtime_error("selected executor does not support domain initialization");
		}
		d->executors.emplace_back(new BackendExecutor(domain, backend, op->GetTimestep()));
		executors.push_back(d->executors.back().get());
		std::cout << "Domain " << domain.id << ": owns X[" << domain.begin << "," << domain.end
		          << "), storage nodes " << domain.StorageNodes() << ", " << backend->GetName()
		          << std::endl;
	}
	d->lower_boundary.resize(shape.Plane());
	d->coordinator.reset(new FDTDCoordinator(*d->plan, d->transport, d->collectives, executors));
	d->coordinator->Bootstrap();
	d->MirrorLowerCurrent();
	return true;
}

void GPU_Backend_Domains::UpdateVoltages()
{
	d->coordinator->BeginPhase(FieldPhase::Voltage);
}
void GPU_Backend_Domains::UpdateCurrents()
{
	d->coordinator->BeginPhase(FieldPhase::Current);
}
void GPU_Backend_Domains::FinishVoltagePhase()
{
	d->coordinator->FinishPhase(FieldPhase::Voltage);
}
void GPU_Backend_Domains::FinishCurrentPhase()
{
	d->coordinator->FinishPhase(FieldPhase::Current);
	d->MirrorLowerCurrent();
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
	std::unique_ptr<DomainExtension> result(new DomainExtension());
	for (auto& backend : d->devices)
	{
		backend->Activate();
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
#endif

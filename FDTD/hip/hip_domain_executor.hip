// SPDX-License-Identifier: GPL-3.0-or-later
// Cartesian unfused curls follow GPU_Backend_Reference / Engine at 370da07.
#include "hip_domain_executor.h"
#include <hip/hip_runtime.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>

namespace openems
{
namespace distributed
{
namespace
{
void Check(hipError_t e, const char* operation)
{
	if (e != hipSuccess)
	{
		throw std::runtime_error(std::string(operation) + ": " + hipGetErrorString(e));
	}
}

struct DeviceGuard
{
	int old;

	explicit DeviceGuard(int device)
	{
		Check(hipGetDevice(&old), "get device");
		Check(hipSetDevice(device), "select domain device");
	}

	~DeviceGuard()
	{
		hipSetDevice(old);
	}
};

struct KernelDomain
{
	uint32_t nx, ny, nz, begin, count;
	size_t plane, storage, owned;
};

__global__ void VoltageKernel(
    float* v, const float* i, const float* vv, const float* vi, KernelDomain d)
{
	const size_t j = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
	if (j >= d.owned)
	{
		return;
	}
	const uint32_t x = uint32_t(j / d.plane) + 1, y = uint32_t(j / d.nz) % d.ny,
	               z = uint32_t(j % d.nz);
	const size_t p = d.plane, s = d.storage, k = j + p, sx = (d.begin + x - 1) ? p : 0,
	             sy = y ? d.nz : 0, sz = z ? 1 : 0;
	v[k] *= vv[j];
	v[k] += vi[j] * (i[2 * s + k] - i[2 * s + k - sy] - i[s + k] + i[s + k - sz]);
	v[s + k] *= vv[d.owned + j];
	v[s + k] += vi[d.owned + j] * (i[k] - i[k - sz] - i[2 * s + k] + i[2 * s + k - sx]);
	v[2 * s + k] *= vv[2 * d.owned + j];
	v[2 * s + k] += vi[2 * d.owned + j] * (i[s + k] - i[s + k - sx] - i[k] + i[k - sy]);
}

__global__ void CurrentKernel(
    const float* v, float* i, const float* ii, const float* iv, KernelDomain d)
{
	const size_t j = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
	if (j >= d.owned)
	{
		return;
	}
	const uint32_t x = uint32_t(j / d.plane) + 1, y = uint32_t(j / d.nz) % d.ny,
	               z = uint32_t(j % d.nz);
	if (d.begin + x - 1 == d.nx - 1 || y == d.ny - 1 || z == d.nz - 1)
	{
		return;
	}
	const size_t p = d.plane, s = d.storage, k = j + p;
	i[k] *= ii[j];
	i[k] += iv[j] * (v[2 * s + k] - v[2 * s + k + d.nz] - v[s + k] + v[s + k + 1]);
	i[s + k] *= ii[d.owned + j];
	i[s + k] += iv[d.owned + j] * (v[k] - v[k + 1] - v[2 * s + k] + v[2 * s + k + p]);
	i[2 * s + k] *= ii[2 * d.owned + j];
	i[2 * s + k] += iv[2 * d.owned + j] * (v[s + k] - v[s + k + p] - v[k] + v[k + d.nz]);
}

struct Pinned
{
	float* data = nullptr;

	explicit Pinned(size_t count)
	{
		Check(hipHostMalloc(reinterpret_cast<void**>(&data),
		          CheckedMultiply(count, sizeof(float)),
		          hipHostMallocPortable),
		    "allocate pinned halo");
	}

	~Pinned()
	{
		if (data)
		{
			hipHostFree(data);
		}
	}
};

class EventCompletion : public Completion
{
public:
	EventCompletion(int device,
	    std::shared_ptr<void> owner,
	    hipStream_t stream,
	    std::function<void()> finish = {})
	    : m_Device(device), m_Owner(owner), m_Event(nullptr), m_Finish(finish), m_Done(false)
	{
		DeviceGuard guard(device);
		Check(hipEventCreateWithFlags(&m_Event, hipEventDisableTiming), "create domain event");
		try
		{
			Check(hipEventRecord(m_Event, stream), "record domain event");
		}
		catch (...)
		{
			hipEventDestroy(m_Event);
			throw;
		}
	}

	~EventCompletion()
	{
		if (!m_Event)
		{
			return;
		}
		try
		{
			Wait();
			DeviceGuard guard(m_Device);
			hipEventDestroy(m_Event);
		}
		catch (const std::exception& e)
		{
			std::cerr << "HIP completion cleanup: " << e.what() << '\n';
		}
	}

	bool Ready() override
	{
		if (m_Done)
		{
			return true;
		}
		DeviceGuard guard(m_Device);
		auto status = hipEventQuery(m_Event);
		if (status == hipErrorNotReady)
		{
			return false;
		}
		Check(status, "query domain event");
		Finish();
		return true;
	}

	void Wait() override
	{
		if (m_Done)
		{
			return;
		}
		DeviceGuard guard(m_Device);
		Check(hipEventSynchronize(m_Event), "wait for domain event");
		Finish();
	}

private:
	void Finish()
	{
		if (m_Finish)
		{
			m_Finish();
		}
		m_Finish = {};
		m_Done = true;
	}

	int m_Device;
	std::shared_ptr<void> m_Owner;
	hipEvent_t m_Event;
	std::function<void()> m_Finish;
	bool m_Done;
};
}

struct HIPDomainExecutor::Impl
{
	Domain domain;
	double timestep;
	int device;
	std::string name;
	KernelDomain kernel;
	hipStream_t stream = nullptr;
	std::array<float*, 6> arrays = { { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr } };
	std::map<uint32_t, std::shared_ptr<Pinned>> send, receive;
	size_t bytes = 0;

	Impl(const LocalOperatorData& data, int selected)
	    : domain(data.domain), timestep(data.timestep), device(selected)
	{
		data.Validate();
		DeviceGuard guard(device);
		hipDeviceProp_t prop;
		Check(hipGetDeviceProperties(&prop, device), "get selected GPU properties");
		name = prop.name;
		kernel = { domain.grid.nx,
			domain.grid.ny,
			domain.grid.nz,
			domain.begin,
			domain.end - domain.begin,
			domain.grid.Plane(),
			domain.StorageNodes(),
			domain.OwnedNodes() };
		if ((kernel.owned + 255) / 256 > uint64_t(prop.maxGridSize[0]))
		{
			throw std::invalid_argument("domain exceeds launch range");
		}
		size_t free = 0, total = 0;
		Check(hipMemGetInfo(&free, &total), "query domain device memory");
		bytes = BasicDomainBytes(domain);
		if (free < CheckedAdd(bytes, 256 * 1024 * 1024))
		{
			throw std::runtime_error("insufficient device memory with reserve");
		}
		try
		{
			Check(hipStreamCreateWithFlags(&stream, hipStreamNonBlocking), "create domain stream");
			for (size_t n = 0; n < 6; ++n)
			{
				const size_t count = CheckedMultiply(3, n < 2 ? kernel.storage : kernel.owned);
				Check(hipMalloc(reinterpret_cast<void**>(&arrays[n]),
				          CheckedMultiply(count, sizeof(float))),
				    "allocate local domain array");
				if (n < 2)
				{
					Check(hipMemsetAsync(arrays[n], 0, count * sizeof(float), stream),
					    "zero local fields");
				}
				else
				{
					Check(hipMemcpyAsync(arrays[n],
					          data.coefficients[n - 2].data(),
					          count * sizeof(float),
					          hipMemcpyHostToDevice,
					          stream),
					    "upload coefficient slice");
				}
			}
			for (int f = 0; f < 2; ++f)
			{
				if (!data.initial_fields[f].empty())
				{
					for (uint32_t c = 0; c < 3; ++c)
					{
						Check(hipMemcpyAsync(arrays[f] + domain.LocalIndex(c, 1, 0, 0),
						          data.initial_fields[f].data() + c * kernel.owned,
						          kernel.owned * sizeof(float),
						          hipMemcpyHostToDevice,
						          stream),
						    "upload owned initial field");
					}
				}
			}
			Check(hipStreamSynchronize(stream), "initialize domain");
		}
		catch (...)
		{
			Cleanup();
			throw;
		}
	}

	void Cleanup()
	{
		if (stream)
		{
			hipStreamSynchronize(stream);
		}
		for (auto& p : arrays)
		{
			if (p)
			{
				hipFree(p);
			}
			p = nullptr;
		}
		if (stream)
		{
			hipStreamDestroy(stream);
		}
		stream = nullptr;
	}

	~Impl()
	{
		int previous = 0;
		const bool have_previous = hipGetDevice(&previous) == hipSuccess;
		if (hipSetDevice(device) == hipSuccess)
		{
			Cleanup();
		}
		else
		{
			std::cerr << "HIP domain cleanup: device unavailable\n";
		}
		if (have_previous)
		{
			hipSetDevice(previous);
		}
	}
};

HIPDomainExecutor::HIPDomainExecutor(const LocalOperatorData& data, int device)
    : m_Impl(std::make_shared<Impl>(data, device))
{
}

HIPDomainExecutor::~HIPDomainExecutor()
{
}

uint32_t HIPDomainExecutor::Id() const
{
	return m_Impl->domain.id;
}

double HIPDomainExecutor::Timestep() const
{
	return m_Impl->timestep;
}

std::string HIPDomainExecutor::DeviceName() const
{
	return m_Impl->name;
}

size_t HIPDomainExecutor::AllocatedDeviceBytes() const
{
	return m_Impl->bytes;
}

CompletionPtr HIPDomainExecutor::Update(FieldPhase phase, uint64_t)
{
	auto d = m_Impl;
	DeviceGuard guard(d->device);
	const uint32_t blocks = uint32_t((d->kernel.owned + 255) / 256);
	if (phase == FieldPhase::Voltage)
	{
		hipLaunchKernelGGL(VoltageKernel,
		    dim3(blocks),
		    dim3(256),
		    0,
		    d->stream,
		    d->arrays[0],
		    d->arrays[1],
		    d->arrays[2],
		    d->arrays[3],
		    d->kernel);
	}
	else if (phase == FieldPhase::Current)
	{
		hipLaunchKernelGGL(CurrentKernel,
		    dim3(blocks),
		    dim3(256),
		    0,
		    d->stream,
		    d->arrays[0],
		    d->arrays[1],
		    d->arrays[4],
		    d->arrays[5],
		    d->kernel);
	}
	else
	{
		throw std::invalid_argument("invalid field phase");
	}
	Check(hipGetLastError(), "launch domain half-step");
	return std::make_shared<EventCompletion>(d->device, d, d->stream);
}

CompletionPtr HIPDomainExecutor::Pack(const HaloRecipe& r, HostBuffer buffer)
{
	auto d = m_Impl;
	DeviceGuard guard(d->device);
	if (!buffer || buffer->size() != r.Samples() || r.source != Id() ||
	    r.plane_samples != d->domain.grid.Plane() || r.source_x == 0 ||
	    r.source_x > d->domain.end - d->domain.begin)
	{
		throw std::invalid_argument("invalid HIP packing recipe");
	}
	auto& pin = d->send[r.channel];
	if (!pin)
	{
		pin = std::make_shared<Pinned>(r.Samples());
	}
	const auto stage = pin;
	const auto field = d->arrays[r.phase == FieldPhase::Voltage ? 0 : 1];
	for (uint32_t c = 1; c < 3; ++c)
	{
		Check(hipMemcpyAsync(stage->data + (c - 1) * r.plane_samples,
		          field + d->domain.LocalIndex(c, r.source_x, 0, 0),
		          r.plane_samples * sizeof(float),
		          hipMemcpyDeviceToHost,
		          d->stream),
		    "pack domain halo");
	}
	return std::make_shared<EventCompletion>(d->device,
	    d,
	    d->stream,
	    [stage, buffer]()
	{
		std::copy_n(stage->data, buffer->size(), buffer->data());
	});
}

CompletionPtr HIPDomainExecutor::Unpack(const HaloRecipe& r, HostBuffer buffer)
{
	auto d = m_Impl;
	DeviceGuard guard(d->device);
	if (!buffer || buffer->size() != r.Samples() || r.destination != Id() ||
	    r.plane_samples != d->domain.grid.Plane() ||
	    (r.destination_x != 0 && r.destination_x != d->domain.end - d->domain.begin + 1))
	{
		throw std::invalid_argument("invalid HIP unpacking recipe");
	}
	auto& pin = d->receive[r.channel];
	if (!pin)
	{
		pin = std::make_shared<Pinned>(r.Samples());
	}
	const auto stage = pin;
	std::copy(buffer->begin(), buffer->end(), stage->data);
	const auto field = d->arrays[r.phase == FieldPhase::Voltage ? 0 : 1];
	for (uint32_t c = 1; c < 3; ++c)
	{
		Check(hipMemcpyAsync(field + d->domain.LocalIndex(c, r.destination_x, 0, 0),
		          stage->data + (c - 1) * r.plane_samples,
		          r.plane_samples * sizeof(float),
		          hipMemcpyHostToDevice,
		          d->stream),
		    "unpack domain halo");
	}
	return std::make_shared<EventCompletion>(d->device,
	    d,
	    d->stream,
	    [stage]()
	{
	});
}

std::vector<float> HIPDomainExecutor::OwnedField(FieldPhase phase) const
{
	auto d = m_Impl;
	DeviceGuard guard(d->device);
	std::vector<float> result(3 * d->domain.OwnedNodes());
	const auto field = d->arrays[phase == FieldPhase::Voltage ? 0 : 1];
	for (uint32_t c = 0; c < 3; ++c)
	{
		Check(hipMemcpyAsync(result.data() + c * d->domain.OwnedNodes(),
		          field + d->domain.LocalIndex(c, 1, 0, 0),
		          d->domain.OwnedNodes() * sizeof(float),
		          hipMemcpyDeviceToHost,
		          d->stream),
		    "gather owned domain fields");
	}
	Check(hipStreamSynchronize(d->stream), "complete owned field gather");
	return result;
}

}
}

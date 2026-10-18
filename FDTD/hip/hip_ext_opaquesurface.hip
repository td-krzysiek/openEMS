// SPDX-License-Identifier: GPL-3.0-or-later
#include "hip/hip_runtime.h"
#include "hip_internal.h"
#include "FDTD/gpu_opaquesurface.h"
#include "FDTD/extensions/engine_ext_opaquesurface.h"
#include <unordered_set>

struct OpaqueData
{
	const uint64_t* index;
	const uint32_t* offset;
	const float *vv, *vi, *a, *b, *memory;
	float *voltage, *state;
};

__global__ void opaque_apply(float* volt, const float* curr, OpaqueData D, uint32_t count)
{
	const size_t i = size_t(blockIdx.x)*blockDim.x+threadIdx.x;
	if (i>=count) return;
	const uint64_t* ix=D.index+5*i;
	const float curl=curr[ix[1]]-curr[ix[2]]-curr[ix[3]]+curr[ix[4]];
	float mem=0;
	for (uint32_t j=D.offset[i]; j<D.offset[i+1]; ++j)
		mem+=D.memory[j]*D.state[j];
	const float old=D.voltage[i];
	const float next=D.vv[i]*old+D.vi[i]*(curl-mem);
	for (uint32_t j=D.offset[i]; j<D.offset[i+1]; ++j)
		D.state[j]=D.a[j]*D.state[j]+D.b[j]*(next+old);
	D.voltage[i]=next;
	volt[ix[0]]=next;
}

class HIP_Ext_OpaqueSurface : public GPU_Extension
{
public:
	HIP_Ext_OpaqueSurface(GPU_Backend_HIP::Impl* impl, Operator_Ext_OpaqueSurface* op)
	    : d(impl), count(static_cast<uint32_t>(op->edges.size()))
	{
		GPU_OpaqueSurfaceData h(*op,d->dim.nx,d->dim.ny,d->dim.nz);
		if (!count) return;
		data.index=d->Alloc<uint64_t>(h.index.size(),h.index.data());
		data.offset=d->Alloc<uint32_t>(h.offset.size(),h.offset.data());
		data.vv=d->Alloc<float>(h.vv.size(),h.vv.data());
		data.vi=d->Alloc<float>(h.vi.size(),h.vi.data());
		data.a=d->Alloc<float>(h.a.size(),h.a.data());
		data.b=d->Alloc<float>(h.b.size(),h.b.data());
		data.memory=d->Alloc<float>(h.memory.size(),h.memory.data());
		data.voltage=d->Alloc<float>(count);
		data.state=d->Alloc<float>(h.a.size());
		// The fused volume step speculates currents before this voltage update.
		// Register every modified edge so it recomputes the affected currents.
		for (size_t i=0; i<count; ++i) modified.insert(h.index[5*i]);
		d->volt_modified.reserve(d->volt_modified.size()+count);
		for (size_t i=0; i<count; ++i) d->volt_modified.push_back(h.index[5*i]);
	}
	virtual ~HIP_Ext_OpaqueSurface()
	{
		// Also called if another extension causes host fallback during Init().
		// Remove one registration per edge, retaining any other extension's copy.
		d->volt_modified.erase(std::remove_if(d->volt_modified.begin(),d->volt_modified.end(),
		    [this](size_t index) {return modified.erase(index)>0;}),d->volt_modified.end());
	}
	virtual void Apply2Voltages()
	{
		if (count) HIP_Launch(d,"opaque_apply",opaque_apply,count,1,1,d->volt,
		                     static_cast<const float*>(d->curr),data,count);
	}
private:
	GPU_Backend_HIP::Impl* d;
	uint32_t count;
	OpaqueData data;
	std::unordered_set<size_t> modified;
};

GPU_Extension* HIP_CreateExt_OpaqueSurface(GPU_Backend_HIP::Impl* d, Engine_Extension* ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_OpaqueSurface*>(ext)) return NULL;
	auto* op=dynamic_cast<Operator_Ext_OpaqueSurface*>(ext->GetOperatorExtension());
	return op ? new HIP_Ext_OpaqueSurface(d,op) : NULL;
}

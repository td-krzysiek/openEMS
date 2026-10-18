// SPDX-License-Identifier: GPL-3.0-or-later
#include "metal_internal.h"
#include "FDTD/gpu_opaquesurface.h"
#include "FDTD/extensions/engine_ext_opaquesurface.h"

static const char* OPAQUE_SOURCE = R"MSL(
kernel void opaque_apply(device float* volt [[buffer(0)]],
                         const device float* curr [[buffer(1)]],
                         const device ulong* index [[buffer(2)]],
                         const device uint* offset [[buffer(3)]],
                         const device float* vv [[buffer(4)]],
                         const device float* vi [[buffer(5)]],
                         const device float* a [[buffer(6)]],
                         const device float* b [[buffer(7)]],
                         const device float* memory [[buffer(8)]],
                         device float* voltage [[buffer(9)]],
                         device float* state [[buffer(10)]],
                         constant uint& count [[buffer(11)]],
                         uint i [[thread_position_in_grid]])
{
	if (i>=count) return;
	const ulong k=ulong(5)*i;
	const float curl=curr[index[k+1]]-curr[index[k+2]]-curr[index[k+3]]+curr[index[k+4]];
	float mem=0;
	for (uint j=offset[i]; j<offset[i+1]; ++j)
		mem+=memory[j]*state[j];
	const float old=voltage[i];
	const float next=vv[i]*old+vi[i]*(curl-mem);
	for (uint j=offset[i]; j<offset[i+1]; ++j)
		state[j]=a[j]*state[j]+b[j]*(next+old);
	voltage[i]=next;
	volt[index[k]]=next;
}
)MSL";

class Metal_Ext_OpaqueSurface : public GPU_Extension
{
public:
	Metal_Ext_OpaqueSurface(GPU_Backend_Metal::Impl* impl, Operator_Ext_OpaqueSurface* op)
	    : d(impl), count(static_cast<uint32_t>(op->edges.size()))
	{
		GPU_OpaqueSurfaceData h(*op,d->dim.nx,d->dim.ny,d->dim.nz);
		if (!count) return;
		index=d->NewBuffer(h.index.size()*sizeof(uint64_t),h.index.data());
		offset=d->NewBuffer(h.offset.size()*sizeof(uint32_t),h.offset.data());
		vv=d->NewBuffer(h.vv.size()*sizeof(float),h.vv.data());
		vi=d->NewBuffer(h.vi.size()*sizeof(float),h.vi.data());
		a=d->NewBuffer(h.a.size()*sizeof(float),h.a.data());
		b=d->NewBuffer(h.b.size()*sizeof(float),h.b.data());
		memory=d->NewBuffer(h.memory.size()*sizeof(float),h.memory.data());
		voltage=d->NewBuffer(count*sizeof(float));
		state=d->NewBuffer(h.a.size()*sizeof(float));
		d->Pipeline(OPAQUE_SOURCE,"opaque_apply");
	}
	virtual void Apply2Voltages()
	{
		if (!count) return;
		id<MTLComputePipelineState> pso=d->Pipeline(OPAQUE_SOURCE,"opaque_apply");
		id<MTLComputeCommandEncoder> enc=d->Encoder();
		[enc setComputePipelineState:pso];
		[enc setBuffer:d->volt offset:0 atIndex:0];
		[enc setBuffer:d->curr offset:0 atIndex:1];
		[enc setBuffer:index offset:0 atIndex:2];
		[enc setBuffer:offset offset:0 atIndex:3];
		[enc setBuffer:vv offset:0 atIndex:4];
		[enc setBuffer:vi offset:0 atIndex:5];
		[enc setBuffer:a offset:0 atIndex:6];
		[enc setBuffer:b offset:0 atIndex:7];
		[enc setBuffer:memory offset:0 atIndex:8];
		[enc setBuffer:voltage offset:0 atIndex:9];
		[enc setBuffer:state offset:0 atIndex:10];
		[enc setBytes:&count length:sizeof(count) atIndex:11];
		d->Dispatch(pso,count);
	}
private:
	GPU_Backend_Metal::Impl* d;
	uint32_t count;
	id<MTLBuffer> index, offset, vv, vi, a, b, memory, voltage, state;
};

GPU_Extension* Metal_CreateExt_OpaqueSurface(GPU_Backend_Metal::Impl* d, Engine_Extension* ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_OpaqueSurface*>(ext)) return NULL;
	auto* op=dynamic_cast<Operator_Ext_OpaqueSurface*>(ext->GetOperatorExtension());
	return op ? new Metal_Ext_OpaqueSurface(d,op) : NULL;
}

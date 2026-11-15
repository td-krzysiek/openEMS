/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "metal_internal.h"
#include "FDTD/engine.h"
#include "FDTD/gpu_dispersive_domain.h"
#include "FDTD/gpu_lossy_fusion.h"
#include "FDTD/extensions/operator_ext_lossymetal.h"
#include <cstdlib>
#include <iostream>
#include "FDTD/extensions/operator_ext_lorentzmaterial.h"
#include "FDTD/extensions/engine_ext_lorentzmaterial.h"

// Drude/Lorentz materials and conducting sheets, see Engine_Ext_LorentzMaterial
// and Engine_Ext_Dispersive. One dispatch per dispersion order, one thread per mesh
// position of that order (unique within an order). The ADE state and coefficients
// are stored per order as [direction][position].
static const char* LORENTZ_SOURCE = R"MSL(
struct ADEParam { uint count; uint lorentz; uint sn; };

// ADE update before the main update, reads the field
kernel void lorentz_pre(const device float* field  [[buffer(0)]],
                        device float* ade          [[buffer(1)]],
                        device float* lor          [[buffer(2)]],
                        const device float* c_int  [[buffer(3)]],
                        const device float* c_ext  [[buffer(4)]],
                        const device float* c_lor  [[buffer(5)]],
                        const device uint* pos     [[buffer(6)]],
                        constant ADEParam& P       [[buffer(7)]],
                        uint i [[thread_position_in_grid]])
{
	if (i>=P.count)
		return;
	for (uint n=0; n<3; ++n)
	{
		const uint k = n*P.count + i;
		const float f = field[n*P.sn + pos[i]];
		if (P.lorentz)
		{
			lor[k] = lor[k] + c_lor[k]*ade[k];
			ade[k] = ade[k] * c_int[k];
			ade[k] = ade[k] + c_ext[k]*(f - lor[k]);
		}
		else
		{
			ade[k] = ade[k] * c_int[k];
			ade[k] = ade[k] + c_ext[k]*f;
		}
	}
}

// subtract the ADE current from the field
kernel void dispersive_apply(device float* field      [[buffer(0)]],
                             const device float* ade  [[buffer(1)]],
                             const device uint* pos   [[buffer(6)]],
                             constant ADEParam& P     [[buffer(7)]],
                             uint i [[thread_position_in_grid]])
{
	if (i>=P.count)
		return;
	for (uint n=0; n<3; ++n)
	{
		const uint g = n*P.sn + pos[i];
		field[g] = field[g] - ade[n*P.count + i];
	}
}

kernel void lossy_pre(const device float* field [[buffer(0)]],
    const device uint* pos [[buffer(1)]],
    constant ADEParam& P [[buffer(2)]],
    device float* ade0 [[buffer(3)]],
    device float* ade1 [[buffer(4)]],
    device float* ade2 [[buffer(5)]],
    device float* ade3 [[buffer(6)]],
    device float* ade4 [[buffer(7)]],
    device float* ade5 [[buffer(8)]],
    device float* ade6 [[buffer(9)]],
    device float* ade7 [[buffer(10)]],
    const device float* c_int0 [[buffer(11)]],
    const device float* c_int1 [[buffer(12)]],
    const device float* c_int2 [[buffer(13)]],
    const device float* c_int3 [[buffer(14)]],
    const device float* c_int4 [[buffer(15)]],
    const device float* c_int5 [[buffer(16)]],
    const device float* c_int6 [[buffer(17)]],
    const device float* c_int7 [[buffer(18)]],
    const device float* c_ext0 [[buffer(19)]],
    const device float* c_ext1 [[buffer(20)]],
    const device float* c_ext2 [[buffer(21)]],
    const device float* c_ext3 [[buffer(22)]],
    const device float* c_ext4 [[buffer(23)]],
    const device float* c_ext5 [[buffer(24)]],
    const device float* c_ext6 [[buffer(25)]],
    const device float* c_ext7 [[buffer(26)]],
    uint i [[thread_position_in_grid]])
{
	if (i >= P.count) return;
	device float* ade[8] = {ade0, ade1, ade2, ade3, ade4, ade5, ade6, ade7};
	const device float* c_int[8] = {c_int0, c_int1, c_int2, c_int3, c_int4, c_int5, c_int6, c_int7};
	const device float* c_ext[8] = {c_ext0, c_ext1, c_ext2, c_ext3, c_ext4, c_ext5, c_ext6, c_ext7};
	for (uint n = 0; n < 3; ++n)
	{
		const uint g = n*P.sn + pos[i];
		const uint k = n*P.count + i;
		const float f = field[g];
		for (uint o = 0; o < 8; ++o)
		{
			ade[o][k] = ade[o][k] * c_int[o][k];
			ade[o][k] = ade[o][k] + c_ext[o][k]*f;
		}
	}
}

kernel void lossy_apply(device float* field [[buffer(0)]],
    const device uint* pos [[buffer(1)]],
    constant ADEParam& P [[buffer(2)]],
    const device float* ade0 [[buffer(3)]],
    const device float* ade1 [[buffer(4)]],
    const device float* ade2 [[buffer(5)]],
    const device float* ade3 [[buffer(6)]],
    const device float* ade4 [[buffer(7)]],
    const device float* ade5 [[buffer(8)]],
    const device float* ade6 [[buffer(9)]],
    const device float* ade7 [[buffer(10)]],
    uint i [[thread_position_in_grid]])
{
	if (i >= P.count) return;
	const device float* ade[8] = {ade0, ade1, ade2, ade3, ade4, ade5, ade6, ade7};
	for (uint n = 0; n < 3; ++n)
	{
		const uint g = n*P.sn + pos[i];
		const uint k = n*P.count + i;
		float f = field[g];
		// Preserve branch-order FP32 subtraction.
		for (uint o = 0; o < 8; ++o) f = f - ade[o][k];
		field[g] = f;
	}
}
)MSL";

class Metal_Ext_LorentzMaterial : public GPU_Extension
{
public:
	Metal_Ext_LorentzMaterial(GPU_Backend_Metal::Impl* impl, Operator_Ext_LorentzMaterial* op_ext);

	virtual void DoPreVoltageUpdates()
	{
		if (m_Fusion == GPU_LossyFusion::All) Fused(true);
		else for (auto& order : m_Volt) Pre(order, d->volt);
	}
	virtual void Apply2Voltages()
	{
		if (m_Fusion != GPU_LossyFusion::Off) Fused(false);
		else for (auto& order : m_Volt) Apply(order, d->volt);
	}
	virtual void DoPreCurrentUpdates() {for (size_t o=0; o<m_Curr.size(); ++o) Pre(m_Curr[o], d->curr);}
	virtual void Apply2Current()       {for (size_t o=0; o<m_Curr.size(); ++o) Apply(m_Curr[o], d->curr);}

protected:
	//! one dispersion order of the voltages or currents
	struct Order
	{
		struct {uint32_t count, lorentz, sn;} param;
		id<MTLBuffer> pos;          //!< flat NIJK index of direction 0
		id<MTLBuffer> ade, lor;     //!< state
		id<MTLBuffer> c_int, c_ext, c_lor;
	};

	void Setup(std::vector<Order>& orders, int o, unsigned int count, unsigned int** pos,
	           bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor);
	void Pre(Order& o, id<MTLBuffer> field);
	void Apply(Order& o, id<MTLBuffer> field);
	void Fused(bool pre);

	GPU_Backend_Metal::Impl* d;
	std::vector<Order> m_Volt;
	std::vector<Order> m_Curr;
	GPU_LossyFusion m_Fusion = GPU_LossyFusion::Off;
};

Metal_Ext_LorentzMaterial::Metal_Ext_LorentzMaterial(GPU_Backend_Metal::Impl* impl, Operator_Ext_LorentzMaterial* op_ext)
{
	d = impl;
	for (int o=0; o<op_ext->m_Order; ++o)
	{
		const unsigned int count = op_ext->m_LM_Count.at(o);
		if (op_ext->m_volt_ADE_On[o])
			Setup(m_Volt, o, count, op_ext->m_LM_pos[o], op_ext->m_volt_Lor_ADE_On[o],
			      op_ext->v_int_ADE[o], op_ext->v_ext_ADE[o], op_ext->m_volt_Lor_ADE_On[o] ? op_ext->v_Lor_ADE[o] : NULL);
		if (op_ext->m_curr_ADE_On[o])
			Setup(m_Curr, o, count, op_ext->m_LM_pos[o], op_ext->m_curr_Lor_ADE_On[o],
			      op_ext->i_int_ADE[o], op_ext->i_ext_ADE[o], op_ext->m_curr_Lor_ADE_On[o] ? op_ext->i_Lor_ADE[o] : NULL);
	}
	d->Pipeline(LORENTZ_SOURCE, "lorentz_pre");
	d->Pipeline(LORENTZ_SOURCE, "dispersive_apply");
	const auto requested = GPU_ParseLossyFusion(std::getenv("OPENEMS_GPU_LOSSY_FUSION"));
	if (requested != GPU_LossyFusion::Off)
	{
		std::vector<GPU_LossyBranch> branches;
		for (int o = 0; o < op_ext->m_Order; ++o)
			branches.push_back({op_ext->m_LM_Count.at(o), op_ext->m_LM_pos[o],
			    op_ext->m_volt_ADE_On[o], op_ext->m_curr_ADE_On[o], op_ext->m_volt_Lor_ADE_On[o]});
		if (GPU_CanFuseLossy(dynamic_cast<Operator_Ext_LossyMetal*>(op_ext) != nullptr, branches) &&
		    m_Volt.size() == 8 && m_Curr.empty())
		{
			m_Fusion = requested;
			d->Pipeline(LORENTZ_SOURCE, "lossy_apply");
			if (requested == GPU_LossyFusion::All) d->Pipeline(LORENTZ_SOURCE, "lossy_pre");
			if (g_settings.GetVerboseLevel() > 1)
				std::cout << "GPU LossyMetal branch fusion: " << (requested == GPU_LossyFusion::All ? "all" : "apply")
				          << ", owned positions " << m_Volt[0].param.count << std::endl;
		}
	}
}


void Metal_Ext_LorentzMaterial::Setup(std::vector<Order>& orders, int o, unsigned int count, unsigned int** pos,
                                      bool lorentz, FDTD_FLOAT** c_int, FDTD_FLOAT** c_ext, FDTD_FLOAT** c_lor)
{
	UNUSED(o);
	GPU_DispersiveDomain local(d->domain, count, pos, d->dim.ny, d->dim.nz);
	count = static_cast<unsigned int>(local.source.size());
	if (count == 0)
	{
		return;
	}
	Order ord;
	ord.param.count = count;
	ord.param.lorentz = lorentz;
	ord.param.sn = d->numCells;

	const std::vector<uint32_t> flat(local.position.begin(), local.position.end());
	ord.pos = d->NewBuffer(count*sizeof(uint32_t), flat.data());

	ord.ade = d->NewBuffer(3*(size_t)count*sizeof(float));
	ord.lor = d->NewBuffer(3*(size_t)count*sizeof(float));
	auto upload = [&](FDTD_FLOAT** values)
	{
		const auto data = local.Coefficients(values);
		return d->NewBuffer(data.size()*sizeof(float), data.data());
	};
	ord.c_int = upload(c_int);
	ord.c_ext = upload(c_ext);
	ord.c_lor = upload(c_lor);
	orders.push_back(ord);
}

void Metal_Ext_LorentzMaterial::Pre(Order& o, id<MTLBuffer> field)
{
	id<MTLComputePipelineState> pso = d->Pipeline(LORENTZ_SOURCE, "lorentz_pre");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:o.ade offset:0 atIndex:1];
	[enc setBuffer:o.lor offset:0 atIndex:2];
	[enc setBuffer:o.c_int offset:0 atIndex:3];
	[enc setBuffer:o.c_ext offset:0 atIndex:4];
	[enc setBuffer:o.c_lor offset:0 atIndex:5];
	[enc setBuffer:o.pos offset:0 atIndex:6];
	[enc setBytes:&o.param length:sizeof(o.param) atIndex:7];
	d->Dispatch(pso, o.param.count);
}

void Metal_Ext_LorentzMaterial::Apply(Order& o, id<MTLBuffer> field)
{
	id<MTLComputePipelineState> pso = d->Pipeline(LORENTZ_SOURCE, "dispersive_apply");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:field offset:0 atIndex:0];
	[enc setBuffer:o.ade offset:0 atIndex:1];
	[enc setBuffer:o.pos offset:0 atIndex:6];
	[enc setBytes:&o.param length:sizeof(o.param) atIndex:7];
	d->Dispatch(pso, o.param.count);
}

void Metal_Ext_LorentzMaterial::Fused(bool pre)
{
	id<MTLComputePipelineState> pso = d->Pipeline(LORENTZ_SOURCE, pre ? "lossy_pre" : "lossy_apply");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:m_Volt[0].pos offset:0 atIndex:1];
	[enc setBytes:&m_Volt[0].param length:sizeof(m_Volt[0].param) atIndex:2];
	for (size_t o = 0; o < 8; ++o)
	{
		[enc setBuffer:m_Volt[o].ade offset:0 atIndex:3+o];
		if (pre)
		{
			[enc setBuffer:m_Volt[o].c_int offset:0 atIndex:11+o];
			[enc setBuffer:m_Volt[o].c_ext offset:0 atIndex:19+o];
		}
	}
	d->Dispatch(pso, m_Volt[0].param.count);
}

GPU_Extension* Metal_CreateExt_LorentzMaterial(GPU_Backend_Metal::Impl* d, Engine_Extension* eng_ext, Engine* eng)
{
	UNUSED(eng);
	if (!dynamic_cast<Engine_Ext_LorentzMaterial*>(eng_ext))
		return NULL;
	Operator_Ext_LorentzMaterial* op_ext = dynamic_cast<Operator_Ext_LorentzMaterial*>(eng_ext->GetOperatorExtension());
	if (!op_ext)
		return NULL;
	return new Metal_Ext_LorentzMaterial(d, op_ext);
}

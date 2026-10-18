// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef GPU_OPAQUESURFACE_H
#define GPU_OPAQUESURFACE_H

#include "extensions/operator_ext_opaquesurface.h"
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

// Shared packing for Metal and HIP. Each edge owns its voltage history and
// a variable-length, disjoint slice of passive R-L branch states. Device
// arithmetic is FP32, like the volume fields; the CPU boundary retains FP64.
struct GPU_OpaqueSurfaceData
{
	std::vector<uint64_t> index; // voltage, Hq, Hq(low p), Hp, Hp(low q)
	std::vector<uint32_t> offset;
	std::vector<float> vv, vi, a, b, memory;

	static float Coefficient(double v)
	{
		const float f = static_cast<float>(v);
		if (!std::isfinite(f) || (v!=0 && f==0))
			throw std::runtime_error("Opaque surface coefficient is outside GPU FP32 range");
		return f;
	}

	GPU_OpaqueSurfaceData(const Operator_Ext_OpaqueSurface& op,
	                      uint32_t nx, uint32_t ny, uint32_t nz)
	{
		if (op.edges.size()>std::numeric_limits<uint32_t>::max())
			throw std::runtime_error("Too many opaque GPU edges");
		const uint64_t stride[3] = {uint64_t(ny)*nz, nz, 1};
		const uint64_t cells = uint64_t(nx)*ny*nz;
		offset.push_back(0);
		for (const auto& e : op.edges)
		{
			const auto& fit = op.fits.at(e.fit);
			if (a.size()>std::numeric_limits<uint32_t>::max()-fit.w.size())
				throw std::runtime_error("Too many opaque GPU branches");
			const unsigned int p=(e.n+1)%3, q=(e.n+2)%3;
			const uint64_t node=uint64_t(e.pos[0])*stride[0]+uint64_t(e.pos[1])*stride[1]+e.pos[2];
			index.push_back(e.n*cells+node);
			index.push_back(q*cells+node);
			index.push_back(q*cells+node-(e.pos[p] ? stride[p] : 0));
			index.push_back(p*cells+node);
			index.push_back(p*cells+node-(e.pos[q] ? stride[q] : 0));
			vv.push_back(Coefficient(e.vv)); vi.push_back(Coefficient(e.vi));
			for (size_t j=0; j<fit.w.size(); ++j)
			{
				const double aj=(2*fit.tau[j]-op.dt)/(2*fit.tau[j]+op.dt);
				a.push_back(Coefficient(aj));
				b.push_back(Coefficient(op.dt*fit.w[j]*e.scale/(2*fit.tau[j]+op.dt)));
				memory.push_back(Coefficient((1+aj)/2));
			}
			offset.push_back(static_cast<uint32_t>(a.size()));
		}
	}
};
#endif

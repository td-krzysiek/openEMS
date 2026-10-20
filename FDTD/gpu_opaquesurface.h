// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef GPU_OPAQUESURFACE_H
#define GPU_OPAQUESURFACE_H

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
	std::vector<uint32_t> offset, fit_offset; // per-edge state and shared-fit starts
	std::vector<float> vv, vi, a, b, memory;

	static float Coefficient(double v)
	{
		const float f = static_cast<float>(v);
		if (!std::isfinite(f) || (v!=0 && f==0))
			throw std::runtime_error("Opaque surface coefficient is outside GPU FP32 range");
		return f;
	}

	// The operator supplies fits, edges and dt; no native engine is needed to test packing.
	template<class Operator>
	GPU_OpaqueSurfaceData(const Operator& op,
	                      uint32_t nx, uint32_t ny, uint32_t nz)
	{
		if (op.edges.size()>std::numeric_limits<uint32_t>::max())
			throw std::runtime_error("Too many opaque GPU edges");
		const uint64_t stride[3] = {uint64_t(ny)*nz, nz, 1};
		const uint64_t cells = uint64_t(nx)*ny*nz;
		std::vector<uint32_t> fit_start(op.fits.size(),std::numeric_limits<uint32_t>::max());
		offset.push_back(0);
		for (const auto& e : op.edges)
		{
			const auto& fit = op.fits.at(e.fit);
			if (fit.w.empty() || fit.w.size()!=fit.tau.size() ||
			    fit.w.size()>std::numeric_limits<uint32_t>::max() ||
			    b.size()>std::numeric_limits<uint32_t>::max()-fit.w.size())
				throw std::runtime_error("Too many opaque GPU branches");
			if (fit_start[e.fit]==std::numeric_limits<uint32_t>::max())
			{
				fit_start[e.fit]=static_cast<uint32_t>(a.size());
				for (size_t j=0; j<fit.w.size(); ++j)
				{
					const double aj=(2*fit.tau[j]-op.dt)/(2*fit.tau[j]+op.dt);
					a.push_back(Coefficient(aj));
					// Preserve rounding from double, not from the already rounded a[j].
					memory.push_back(Coefficient((1+aj)/2));
				}
			}
			fit_offset.push_back(fit_start[e.fit]);
			const unsigned int p=(e.n+1)%3, q=(e.n+2)%3;
			const uint64_t node=uint64_t(e.pos[0])*stride[0]+uint64_t(e.pos[1])*stride[1]+e.pos[2];
			index.push_back(e.n*cells+node);
			index.push_back(q*cells+node);
			index.push_back(q*cells+node-(e.pos[p] ? stride[p] : 0));
			index.push_back(p*cells+node);
			index.push_back(p*cells+node-(e.pos[q] ? stride[q] : 0));
			vv.push_back(Coefficient(e.vv)); vi.push_back(Coefficient(e.vi));
			for (size_t j=0; j<fit.w.size(); ++j)
				b.push_back(Coefficient(op.dt*fit.w[j]*e.scale/(2*fit.tau[j]+op.dt)));
			offset.push_back(static_cast<uint32_t>(b.size()));
		}
	}
};
#endif

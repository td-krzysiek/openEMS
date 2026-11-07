// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef GPU_OPAQUESURFACE_H
#define GPU_OPAQUESURFACE_H

#include <cmath>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

// Shared packing for Metal and HIP. Each edge owns its voltage history and
// a variable-length, disjoint slice of passive R-L branch states. Device
// arithmetic is FP32, like the volume fields; the CPU boundary retains FP64.
struct GPU_OpaqueSurfaceData
{
	std::vector<uint64_t> index; // voltage, Hq, Hq(low p), Hp, Hp(low q)
	std::vector<uint32_t> offset, fit_offset; // per-edge state and shared-fit starts
	std::vector<uint32_t> b_index; // optional shared-vector ID (or mixed-order start)
	std::vector<float> vv, vi, a, b, memory;
	// Uniform-order histories are pole-major within compact 32-edge tiles.
	// A partial last tile has no padding. Mixed orders retain edge-major storage.
	enum { TileWidth = 32 };
	uint32_t pole_count = 0;
	uint32_t coefficient_count = 0; // zero retains the original per-edge b layout

	uint32_t CoefficientIndex(uint32_t edge, uint32_t pole) const
	{
		if (!coefficient_count) return StateIndex(edge, pole);
		if (!pole_count) return b_index[edge] + pole;
		const uint32_t id = b_index[edge], first = id - id % TileWidth;
		const uint32_t lanes = std::min(uint32_t(TileWidth), coefficient_count - first);
		return first * pole_count + pole * lanes + id - first;
	}

	uint32_t StateIndex(uint32_t edge, uint32_t pole) const
	{
		if (!pole_count)
		{
			return offset[edge] + pole;
		}
		const uint32_t first = edge - edge % TileWidth;
		const uint32_t lanes = std::min(uint32_t(TileWidth), uint32_t(vv.size()) - first);
		return first * pole_count + pole * lanes + edge - first;
	}

	static float Coefficient(double v)
	{
		const float f = static_cast<float>(v);
		if (!std::isfinite(f) || (v!=0 && f==0))
		{
			throw std::runtime_error("Opaque surface coefficient is outside GPU FP32 range");
		}
		return f;
	}

	// The operator supplies fits, edges and dt; no native engine is needed to test packing.
	template<class Operator>
	GPU_OpaqueSurfaceData(const Operator& op,
	                      uint32_t nx, uint32_t ny, uint32_t nz)
	{
		if (op.edges.size()>std::numeric_limits<uint32_t>::max())
		{
			throw std::runtime_error("Too many opaque GPU edges");
		}
		const uint64_t stride[3] = {uint64_t(ny)*nz, nz, 1};
		const uint64_t cells = uint64_t(nx)*ny*nz;
		std::vector<uint32_t> fit_start(op.fits.size(),std::numeric_limits<uint32_t>::max());
		offset.push_back(0);
		bool uniform = true;
		for (const auto& e : op.edges)
		{
			const auto& fit = op.fits.at(e.fit);
			if (fit.w.empty() || fit.w.size()!=fit.tau.size() ||
			    fit.w.size()>std::numeric_limits<uint32_t>::max() ||
			    offset.back()>std::numeric_limits<uint32_t>::max()-fit.w.size())
			{
				throw std::runtime_error("Too many opaque GPU branches");
			}
			if (offset.size() == 1)
			{
				pole_count = static_cast<uint32_t>(fit.w.size());
			}
			uniform = uniform && fit.w.size() == pole_count;
			offset.push_back(offset.back() + static_cast<uint32_t>(fit.w.size()));
		}
		if (!uniform)
		{
			pole_count = 0;
		}
		// Only constants may alias. Every edge retains its original history slice.
		struct Key
		{
			uint64_t scale;
			uint32_t fit;
			bool operator==(const Key& other) const
			{ return scale == other.scale && fit == other.fit; }
		};
		struct Hash
		{
			size_t operator()(const Key& key) const
			{ return std::hash<uint64_t>()(key.scale) ^ (size_t(key.fit) * 0x9e3779b9u); }
		};
		std::vector<uint32_t> representatives, starts(1, 0);
		{
			std::unordered_map<Key, uint32_t, Hash> vectors;
			b_index.reserve(op.edges.size());
			for (uint32_t i = 0; i < op.edges.size(); ++i)
			{
				const auto& e = op.edges[i];
				uint64_t scale;
				static_assert(sizeof(scale) == sizeof(e.scale), "Opaque scale must be binary64");
				std::memcpy(&scale, &e.scale, sizeof(scale));
				const auto entry = vectors.emplace(Key{scale, e.fit}, uint32_t(representatives.size()));
				if (entry.second)
				{
					representatives.push_back(i);
					starts.push_back(starts.back() + uint32_t(op.fits[e.fit].w.size()));
				}
				b_index.push_back(pole_count ? entry.first->second : starts[entry.first->second]);
			}
		}
		// Avoid extra device storage/indirection when sharing cannot save bytes.
		if (uint64_t(starts.back()) + b_index.size() < offset.back())
		{
			coefficient_count = uint32_t(representatives.size());
			b.resize(starts.back());
		}
		else
		{
			std::vector<uint32_t>().swap(b_index);
			b.resize(offset.back());
		}
		vv.reserve(op.edges.size());
		vi.reserve(op.edges.size());
		for (const auto& e : op.edges)
		{
			const auto& fit = op.fits[e.fit];
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
			vv.push_back(Coefficient(e.vv));
			vi.push_back(Coefficient(e.vi));
		}
		const size_t packed = coefficient_count ? representatives.size() : op.edges.size();
		for (uint32_t k = 0; k < packed; ++k)
		{
			const uint32_t i = coefficient_count ? representatives[k] : k;
			const auto& e = op.edges[i];
			const auto& fit = op.fits[e.fit];
			for (size_t j=0; j<fit.w.size(); ++j)
			{
				b[CoefficientIndex(i, static_cast<uint32_t>(j))] =
				    Coefficient(op.dt*fit.w[j]*e.scale/(2*fit.tau[j]+op.dt));
			}
		}
	}
};
#endif

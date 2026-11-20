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

#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "gpu_coeff_sets.h"
#include "gpu_backend.h"
#include "operator.h"
#include "tools/useful.h"

namespace
{
//! bit pattern of the values of a set, see GPU_FindSets()
struct SetKey
{
	uint32_t v[GPU_MAX_SET_WIDTH];
	unsigned int width;
	bool operator==(const SetKey& o) const {return std::memcmp(v, o.v, width*sizeof(uint32_t))==0;}
};

struct SetKeyHash
{
	size_t operator()(const SetKey& s) const
	{
		uint64_t h = 1469598103934665603ULL;   // FNV-1a
		for (unsigned int k=0; k<s.width; ++k)
			h = (h ^ s.v[k]) * 1099511628211ULL;
		return h;
	}
};
}

bool GPU_UseRawMainCoefficients(const char* value)
{
	if (!value || std::strcmp(value, "auto") == 0)
		return false;
	if (std::strcmp(value, "raw") == 0)
		return true;
	throw std::runtime_error("OPENEMS_GPU_MAIN_COEFF must be auto or raw");
}

bool GPU_FindSets(size_t count, unsigned int width, const std::function<void(size_t, float*)>& get, GPU_CoeffSets& sets)
{
	// below 3 values per item a set index (counted as 4 bytes) saves too little
	if ((width<3) || (width>GPU_MAX_SET_WIDTH))
		return false;
	// at most half the size of the full arrays (width floats per item), counting
	// 4 bytes per index: whether 16 bit ones are enough is not known yet
	const size_t max_sets16 = (size_t)std::numeric_limits<uint16_t>::max()+1;
	// and no more sets than a 32 bit index can number
	const size_t max_sets = std::min((2*width*count - 4*count)/(4*width), (size_t)std::numeric_limits<uint32_t>::max());
	sets.index.resize(count);
	sets.table.clear();
	std::unordered_map<SetKey, uint32_t, SetKeyHash> found;

	SetKey prev;
	bool have_prev = false;
	// the values in blocks, each filled on all threads (get() may compute them, Coeffs_EC), then searched in order
	const size_t block = std::min(count, (size_t)1<<20);
	std::vector<float> buffer(block*width);
	const unsigned int numThreads = std::max(1u, std::min(AvailableThreads(), (unsigned int)std::max<size_t>(1, block/4096)));
	for (size_t b0=0; b0<count; b0+=block)
	{
		const size_t b1 = std::min(count, b0+block);
		ParallelFor(b1-b0, numThreads, [&](unsigned int, size_t begin, size_t end) {
			for (size_t i=b0+begin; i<b0+end; ++i)
				get(i, &buffer[(i-b0)*width]);
		});
		for (size_t i=b0; i<b1; ++i)
		{
			SetKey key;
			key.width = width;
			const float* values = &buffer[(i-b0)*width];
			std::memcpy(key.v, values, width*sizeof(float));
			// neighbours mostly share the set
			if (have_prev && (key==prev))
			{
				sets.index[i] = sets.index[i-1];
				continue;
			}
			std::unordered_map<SetKey, uint32_t, SetKeyHash>::const_iterator it = found.find(key);
			if (it==found.end())
			{
				if (found.size()>=max_sets)
					return false;
				it = found.insert(std::make_pair(key, (uint32_t)found.size())).first;
				sets.table.insert(sets.table.end(), values, values+width);
			}
			sets.index[i] = it->second;
			prev = key;
			have_prev = true;
		}
	}
	sets.count = found.size();
	sets.mode = (sets.count<=max_sets16) ? 1 : 2;
	return true;
}

void GPU_NodeCoefficients(const Operator* op, const unsigned int numLines[3], size_t node, float values[12])
{
	const unsigned int x = node/((size_t)numLines[1]*numLines[2]);
	const unsigned int y = (node/numLines[2]) % numLines[1];
	const unsigned int z = node % numLines[2];
	const unsigned int pos[3] = {x, y, z};
	for (int n=0; n<3; ++n)
		for (int f=0; f<2; ++f)
			op->GetCoeffPair(f, n, pos, values[3*(2*f)+n], values[3*(2*f+1)+n]);
}

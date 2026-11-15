// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_GPU_DOMAIN_LAYOUT_H
#define OPENEMS_GPU_DOMAIN_LAYOUT_H
#include "gpu_backend.h"
#include <algorithm>
#include <stdexcept>

// Fixed X ownership with one halo on each side, including physical boundaries.
// No device API or transport belongs in this coordinate conversion.
struct GPU_DomainLayout
{
	bool enabled = false;
	unsigned int begin = 0, end = 0, global[3] = { 0, 0, 0 };

	bool Owns(unsigned int x) const
	{
		return !enabled || (x >= begin && x < end);
	}

	unsigned int LocalX(unsigned int x) const
	{
		return enabled ? x - begin + 1 : x;
	}

	unsigned int LocalLines() const
	{
		return end - begin + 2;
	}

};
#endif

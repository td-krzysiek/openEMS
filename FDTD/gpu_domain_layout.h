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

	void Coefficients(
	    const Operator* op, const unsigned int local[3], size_t node, float* values) const
	{
		if (!enabled)
		{
			GPU_NodeCoefficients(op, local, node, values);
			return;
		}
		const size_t plane = size_t(local[1]) * local[2];
		const unsigned int x = unsigned(node / plane);
		if (x == 0 || x + 1 == local[0])
		{
			std::fill(values, values + 12, 0.f);
			std::fill(values, values + 3, 1.f);
			std::fill(values + 6, values + 9, 1.f);
			return;
		}
		const unsigned int gx = begin + x - 1;
		GPU_NodeCoefficients(op, global, size_t(gx) * plane + node % plane, values);
		// The original current loop excludes the true global last X line.
		if (gx + 1 == global[0])
		{
			std::fill(values + 6, values + 9, 1.f);
			std::fill(values + 9, values + 12, 0.f);
		}
	}
};
#endif

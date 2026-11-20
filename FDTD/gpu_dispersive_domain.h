// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_GPU_DISPERSIVE_DOMAIN_H
#define OPENEMS_GPU_DISPERSIVE_DOMAIN_H

#include "gpu_domain_layout.h"
#include <vector>

// Preserve per-order source indices while assigning histories only to owners.
// Halo fields are exchanged by the domain runtime; they have no local history.
struct GPU_DispersiveDomain
{
	std::vector<size_t> source;
	std::vector<size_t> position;

	GPU_DispersiveDomain(const GPU_DomainLayout& layout, size_t count,
	                     unsigned int** coordinates, unsigned int ny, unsigned int nz)
	{
		for (size_t i = 0; i < count; ++i)
		{
			if (layout.Owns(coordinates[0][i]))
			{
				source.push_back(i);
				position.push_back((size_t(layout.LocalX(coordinates[0][i])) * ny +
				                    coordinates[1][i]) * nz + coordinates[2][i]);
			}
		}
	}

	template <typename T>
	std::vector<float> Coefficients(T** coefficients) const
	{
		std::vector<float> result(3 * source.size(), 0.f);
		if (coefficients)
		{
			for (unsigned int n = 0; n < 3; ++n)
			{
				for (size_t i = 0; i < source.size(); ++i)
				{
					result[n * source.size() + i] = coefficients[n][source[i]];
				}
			}
		}
		return result;
	}
};

#endif

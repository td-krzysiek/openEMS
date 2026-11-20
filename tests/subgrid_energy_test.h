/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_ENERGY_TEST_H
#define SUBGRID_ENERGY_TEST_H
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/gpu_cartesiansubgrid.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

inline void CheckSubgridEnergy(GPU_CartesianSubgrid &device, const CartesianSubgrid &grid,
	unsigned int workers, const char *kind)
{
	const auto before = device.EnergyDownloadBytes();
	const auto energy = device.ReadEnergy();
	const double expected[2] = {grid.ElectricEnergy(), grid.MagneticEnergy()};
	double error = 0;
	for (unsigned int n = 0; n < 2; ++n)
	{
		const double relative = std::abs(energy[n] - expected[n]) / std::max(expected[n], 1e-300);
		if (!std::isfinite(relative) || relative > 2e-6)
		{
			throw std::runtime_error("device field energy differs from complete-state Galerkin oracle");
		}
		error = std::max(error, relative);
	}
	const auto bytes = device.EnergyDownloadBytes() - before;
	if (!bytes || bytes > 24 * workers)
	{
		throw std::runtime_error("device field energy downloaded more than scalar reductions");
	}
	std::cout << "ENERGY_COMPOSITION kind=" << kind << " error=" << error << " downloaded_bytes=" << bytes << '\n';
}
#endif

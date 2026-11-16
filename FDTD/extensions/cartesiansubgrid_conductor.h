/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIANSUBGRID_CONDUCTOR_H
#define CARTESIANSUBGRID_CONDUCTOR_H
#include "operator_ext_cartesiansubgrid_conductivity.h"
#include <functional>

namespace CartesianSubgridConductor
{
//! Positive SI conductivity on an already classified opaque leaf.
using Sampler = std::function<double(std::size_t)>;
struct Prepared
{
	std::vector<Operator_Ext_CartesianSubgrid_Conductivity::SurfaceAdmittance> surfaces;
	std::size_t faces = 0, fits = 0, normalColumns = 0;
	double maximumRealFitError = 0, maximumImaginaryFitError = 0;
};
//! Smooth exposed-face quadrature only, with no shifted reference plane.
Prepared Prepare(const CartesianSubgrid &grid, const Sampler &sample, double maximumFrequency);
}
#endif

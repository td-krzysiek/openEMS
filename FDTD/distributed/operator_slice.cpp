// SPDX-License-Identifier: GPL-3.0-or-later
#include "operator_slice.h"
#include "FDTD/operator.h"
#include "FDTD/operator_cylinder.h"
#include "FDTD/extensions/operator_ext_excitation.h"
#include <stdexcept>

namespace openems
{
namespace distributed
{
LocalOperatorData SliceOperator(const Operator& op, const Domain& domain)
{
	if (dynamic_cast<const Operator_Cylinder*>(&op))
	{
		throw std::invalid_argument("distributed operator adapter supports Cartesian grids only");
	}
	const GridShape shape = {
		op.GetNumberOfLines(0, true), op.GetNumberOfLines(1, true), op.GetNumberOfLines(2, true)
	};
	if (domain.grid.nx != shape.nx || domain.grid.ny != shape.ny || domain.grid.nz != shape.nz)
	{
		throw std::invalid_argument("partition does not describe the built global operator");
	}
	LocalOperatorData data;
	data.domain = domain;
	data.timestep = op.GetTimestep();
	for (size_t n = 0; n < op.GetNumberOfExtentions(); ++n)
	{
		const auto* ext = op.GetExtension(n);
		if (!ext->IsActive())
		{
			continue;
		}
		const auto* excitation = dynamic_cast<const Operator_Ext_Excitation*>(ext);
		if (excitation && excitation->GetVoltCount() == 0 && excitation->GetCurrCount() == 0)
		{
			continue;
		}
		data.required_extensions.push_back(ext->GetExtensionName());
	}
	if (!data.required_extensions.empty())
	{
		throw std::invalid_argument(
		    "unsupported distributed extension: " + data.required_extensions.front());
	}
	for (auto& a : data.coefficients)
	{
		a.resize(CheckedMultiply(3, domain.OwnedNodes()));
	}
	for (uint32_t c = 0; c < 3; ++c)
	{
		for (uint32_t x = domain.begin; x < domain.end; ++x)
		{
			for (uint32_t y = 0; y < shape.ny; ++y)
			{
				for (uint32_t z = 0; z < shape.nz; ++z)
				{
					const size_t i = c * domain.OwnedNodes() +
					    size_t(x - domain.begin) * shape.Plane() + size_t(y) * shape.nz + z;
					data.coefficients[0][i] = op.GetVV(c, x, y, z);
					data.coefficients[1][i] = op.GetVI(c, x, y, z);
					data.coefficients[2][i] = op.GetII(c, x, y, z);
					data.coefficients[3][i] = op.GetIV(c, x, y, z);
				}
			}
		}
	}
	data.Validate();
	return data;
}
}
}

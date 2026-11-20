/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_OBSERVATIONS_H
#define SUBGRID_OBSERVATIONS_H

#include "cartesian_subgrid.h"
#include "gpu_cartesiansubgrid.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

//! Owning setup-only SI E/H samples and existing U/I functionals.
//! Bind once before stepping to a device created from the same unchanged grid.
//! The device/backend must outlive this reader. No complete field download occurs.
class CartesianSubgridObservations
{
  public:
	CartesianSubgridObservations(const CartesianSubgrid &grid,
		const std::vector<CartesianSubgrid::Integral> &functionals)
		: m_Functionals(functionals)
	{
		for (const auto &functional : m_Functionals)
		{
			auto &selection = functional.magnetic ? m_Magnetic : m_Electric;
			const auto count = functional.magnetic ? grid.FaceCount() : grid.EdgeCount();
			for (const auto &term : functional.weights)
			{
				if (term.first >= count || term.first > std::numeric_limits<std::uint32_t>::max() ||
					!std::isfinite(term.second))
				{
					throw std::invalid_argument("CartesianSubgrid: invalid observation functional");
				}
				selection.push_back(static_cast<std::uint32_t>(term.first));
			}
		}
		for (auto *selection : {&m_Electric, &m_Magnetic})
		{
			std::sort(selection->begin(), selection->end());
			selection->erase(std::unique(selection->begin(), selection->end()), selection->end());
		}
		for (auto &functional : m_Functionals)
		{
			const auto &selection = functional.magnetic ? m_Magnetic : m_Electric;
			for (auto &term : functional.weights)
			{
				term.first = std::lower_bound(selection.begin(), selection.end(), term.first) - selection.begin();
			}
		}
	}

	void Bind(GPU_CartesianSubgrid &execution)
	{
		if (m_BindAttempted)
		{
			throw std::logic_error("CartesianSubgrid: observation binding is startup-only and single-use");
		}
		m_BindAttempted = true;
		execution.ConfigureObservations(m_Electric, m_Magnetic);
		m_Execution = &execution;
	}

	void Read(std::vector<double> &values)
	{
		if (!m_Execution)
		{
			throw std::logic_error("CartesianSubgrid: observation device is not bound");
		}
		m_Execution->ReadObservations(m_ElectricValues, m_MagneticValues);
		if (m_ElectricValues.size() != m_Electric.size() || m_MagneticValues.size() != m_Magnetic.size())
		{
			throw std::logic_error("CartesianSubgrid: observation state count changed");
		}
		values.resize(m_Functionals.size());
		for (std::size_t i = 0; i < m_Functionals.size(); ++i)
		{
			const auto &functional = m_Functionals[i];
			const auto &state = functional.magnetic ? m_MagneticValues : m_ElectricValues;
			double value = 0;
			for (const auto &term : functional.weights)
			{
				value += term.second * state.at(term.first);
			}
			if (!std::isfinite(value))
			{
				throw std::runtime_error("CartesianSubgrid: nonfinite physical observation");
			}
			values[i] = value;
		}
	}

	std::size_t SelectedElectricCount() const { return m_Electric.size(); }
	std::size_t SelectedMagneticCount() const { return m_Magnetic.size(); }

  private:
	std::vector<CartesianSubgrid::Integral> m_Functionals;
	std::vector<std::uint32_t> m_Electric, m_Magnetic;
	std::vector<float> m_ElectricValues, m_MagneticValues;
	GPU_CartesianSubgrid *m_Execution = nullptr;
	bool m_BindAttempted = false;
};
#endif

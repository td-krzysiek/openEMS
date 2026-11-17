/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_native_field_gather.h"
#include "engine_interface_fdtd.h"
#include "operator_cartesiansubgrid.h"
#include "engine.h"
#include <cmath>
#include <memory>

namespace
{
const Engine_Interface_FDTD *ValidateGather(const Engine_Interface_FDTD *interface,
	const Operator_CartesianSubgrid *op, const unsigned int count[3], unsigned int *const positions[3])
{
	if (!interface || !op || interface->GetOperator() != op || !count || !positions)
	{
		throw std::invalid_argument("CartesianSubgrid: native gather requires its matching interface/operator");
	}
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (!positions[n] || !count[n] || count[n] > op->GetNumberOfLines(n, true))
		{
			throw std::invalid_argument("CartesianSubgrid: invalid native gather positions/count");
		}
		for (unsigned int i = 0; i < count[n]; ++i)
		{
			if (positions[n][i] >= op->GetNumberOfLines(n, true))
			{
				throw std::invalid_argument("CartesianSubgrid: native gather position is outside the mesh");
			}
		}
	}
	return interface;
}
}

Field_Gather_Subgrid::Field_Gather_Subgrid(const Engine_Interface_FDTD *interface,
	const Operator_CartesianSubgrid *op, bool magnetic, const unsigned int count[3],
	unsigned int *const positions[3], std::size_t budget)
	: Field_Gather_FDTD(ValidateGather(interface, op, count, positions), op, nullptr, nullptr, magnetic, count, positions), m_Local(op), m_EngineIdentity(interface->GetFDTDEngine())
{
	const auto charge = [&](std::size_t number, std::size_t unit)
	{
		if (number > (budget - m_Estimate) / unit)
		{
			throw BudgetExceeded();
		}
		m_Estimate += number * unit;
	};
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int d = 0; d < 3; ++d)
		{
			charge(m_RegionPos[n][d].capacity() + region.lines[n][d].capacity(), sizeof(unsigned int));
		}
		charge(m_PosLines[n].capacity(), sizeof(unsigned int));
	}
	const std::size_t size = region.Size();
	// Include retained raw values, descriptors, registry/reader objects and preparation slack.
	charge(size, sizeof(Sample) + sizeof(float) + 4 * sizeof(CartesianSubgrid::Integral));
	std::vector<CartesianSubgrid::Integral> functions;
	functions.reserve(size);
	m_Samples.reserve(size);
	// Empty boundary stencils still need a nonnull region-snapshot pointer.
	// Otherwise the existing interpolator selects its absent full host arrays.
	if (size == 0)
	{
		charge(1, sizeof(float));
	}
	m_Raw.resize(size ? size : 1);
	m_Lines = std::size_t(count[0]) * count[1];
	m_Mode = interface->GetInterpolationType();
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (auto i : region.lines[n][0])
		{
			for (auto j : region.lines[n][1])
			{
				for (auto k : region.lines[n][2])
				{
					const unsigned int p[3] = {i, j, k};
					const bool outside = magnetic ?
						p[(n + 1) % 3] + 1 >= op->GetNumberOfLines((n + 1) % 3) ||
						p[(n + 2) % 3] + 1 >= op->GetNumberOfLines((n + 2) % 3) :
						p[n] + 1 >= op->GetNumberOfLines(n);
					CartesianSubgrid::Integral functional;
					if (!outside)
					{
						const auto point = op->Position(n, p, magnetic);
						functional = magnetic ? op->Grid().MagneticSample(n, point) :
							op->Grid().ElectricSample(n, point);
					}
					charge(functional.weights.capacity(), 4 * sizeof(functional.weights[0]));
					functions.push_back(std::move(functional));
					m_Samples.push_back({std::numeric_limits<std::size_t>::max(), op->GetEdgeLength(n, p, magnetic)});
				}
			}
		}
	}
	// No registry mutation until the complete candidate fits its predeclared budget.
	for (std::size_t i = 0; i < size; ++i)
	{
		if (!functions[i].weights.empty())
		{
			m_Samples[i].observation = op->RegisterObservation(functions[i]);
		}
	}
	SetSnapshotRegion(0);
}

void Field_Gather_Subgrid::Evaluate(std::size_t first, std::size_t last,
	ArrayLib::ArrayNIJK<float> &field, const float *source) const
{
	if (source || first > last || last > m_Lines || m_EngIf->GetInterpolationType() != m_Mode)
	{
		throw std::invalid_argument("CartesianSubgrid: invalid live native gather range, source or changed interpolation");
	}
	if (first == last)
	{
		return;
	}
	{
		std::lock_guard<std::mutex> lock(m_Mutex);
		const auto epoch = m_EngIf->GetNumberOfTimesteps();
		if (epoch != m_Epoch)
		{
			for (std::size_t i = 0; i < m_Samples.size(); ++i)
			{
				const auto &sample = m_Samples[i];
				m_Raw[i] = sample.observation == std::numeric_limits<std::size_t>::max() ? 0.f :
					static_cast<float>(m_Local->EvaluateObservation(sample.observation) * sample.length);
				if (!std::isfinite(m_Raw[i]))
				{
					throw std::runtime_error("CartesianSubgrid: nonfinite native gather sample");
				}
			}
			if (m_EngIf->GetNumberOfTimesteps() != epoch)
			{
				throw std::logic_error("CartesianSubgrid: native gather cannot overlap an advance");
			}
			m_Epoch = epoch;
			++m_Readbacks;
		}
	}
	// Existing interpolation divides the rounded raw float by its physical edge length.
	Field_Gather_FDTD::Evaluate(first, last, field, m_Raw.data());
}

void Field_Gather_Subgrid::CheckInterface(const Engine_Interface_FDTD *caller) const
{
	if (caller != m_EngIf || caller->GetFDTDOperator() != m_Local || caller->GetOperator() != m_Local ||
		caller->GetFDTDEngine() != m_EngineIdentity || m_Local->GetEngine() != m_EngineIdentity ||
		caller->GetInterpolationType() != m_Mode)
	{
		throw std::invalid_argument("CartesianSubgrid: DFT gather interface/engine/interpolation changed");
	}
}

int Field_Gather_Subgrid::CreateDeviceDFT(const Engine_Interface_FDTD *caller, unsigned int frequencies) const
{
	CheckInterface(caller);
	if (m_DFTId >= 0)
	{
		throw std::logic_error("CartesianSubgrid: native gather DFT is already registered");
	}
	std::vector<std::size_t> observations;
	observations.reserve(m_Samples.size());
	for (const auto &sample : m_Samples)
	{
		observations.push_back(sample.observation);
	}
	const int id = m_Local->CreateFieldDFT(m_H, observations, frequencies);
	if (id >= 0)
	{
		m_DFTId = id;
		m_DFTFrequencies = frequencies;
	}
	return id;
}

void Field_Gather_Subgrid::AccumulateDeviceDFT(const Engine_Interface_FDTD *caller,
	const std::vector<std::complex<float>> &weights) const
{
	CheckInterface(caller);
	m_Local->AccumulateFieldDFT(m_DFTId, weights);
}

bool Field_Gather_Subgrid::DFTOutputMatches(const Engine_Interface_FDTD *caller,
	const ArrayLib::ArrayNIJK<std::complex<float>> &field) const
{
	CheckInterface(caller);
	return field.data() && field.extent(0) == 3 && field.extent(1) == m_PosLines[0].size() &&
		field.extent(2) == m_PosLines[1].size() && field.extent(3) == m_PosLines[2].size();
}

void Field_Gather_Subgrid::ReadDeviceDFT(const Engine_Interface_FDTD *caller, unsigned int frequency,
	ArrayLib::ArrayNIJK<std::complex<float>> &field) const
{
	if (!DFTOutputMatches(caller, field))
	{
		throw std::invalid_argument("CartesianSubgrid: native DFT output shape differs");
	}
	std::vector<std::complex<double>> sums;
	m_Local->ReadFieldDFT(m_DFTId, frequency, sums);
	if (sums.size() != m_Samples.size())
	{
		throw std::logic_error("CartesianSubgrid: reconstructed DFT raw region size differs");
	}
	std::vector<std::complex<float>> raw(sums.empty() ? 1 : sums.size());
	for (std::size_t i = 0; i < sums.size(); ++i)
	{
		raw[i] = std::complex<float>(sums[i] * m_Samples[i].length);
		if (!std::isfinite(raw[i].real()) || !std::isfinite(raw[i].imag()))
		{
			m_Local->InvalidateObservations();
			throw std::runtime_error("CartesianSubgrid: native DFT raw field exceeds float range");
		}
	}
	// Reconstruction now occurs after accumulation. Per-sample raw-float rounding
	// is not reproduced: qualify this bounded derivative, never claim bit identity.
	Field_Gather_FDTD::EvaluateSums(raw.data(), field);
	for (std::size_t i = 0; i < field.size(); ++i)
	{
		if (!std::isfinite(field.data(i).real()) || !std::isfinite(field.data(i).imag()))
		{
			m_Local->InvalidateObservations();
			throw std::runtime_error("CartesianSubgrid: nonfinite interpolated DFT output");
		}
	}
}

std::uint64_t Field_Gather_Subgrid::Readbacks() const
{
	std::lock_guard<std::mutex> lock(m_Mutex);
	return m_Readbacks;
}

Engine_Field_Gather *Operator_CartesianSubgrid::CreateFieldGather(const Engine_Interface_FDTD *interface,
	bool magnetic, const unsigned int count[3], unsigned int *const positions[3]) const
{
	if (!m_Grid || !m_Engine || m_ObservationPlanSealed || m_Engine->GetNumberOfTimesteps() != 0 || !m_FieldGatherBudget)
	{
		return nullptr;
	}
	if (!count || !positions)
	{
		throw std::invalid_argument("CartesianSubgrid: null native gather mesh");
	}
	for (unsigned int n = 0; n < 3; ++n)
	{
		if (!count[n])
		{
			return nullptr;
		}
	}
	try
	{
		auto gather = std::make_unique<Field_Gather_Subgrid>(interface,
			this, magnetic, count, positions, m_FieldGatherBudget);
		m_FieldGatherBudget -= gather->EstimatedStorageBytes();
		m_FieldGatherUsed += gather->EstimatedStorageBytes();
		return gather.release();
	}
	catch (const Field_Gather_Subgrid::BudgetExceeded &)
	{
		return nullptr; // Existing explicit full observation fallback; no partial registration.
	}
}

void Operator_CartesianSubgrid::SetFieldGatherBudget(std::size_t bytes)
{
	if (m_ObservationPlanSealed || m_FieldGatherUsed || (m_Engine && m_Engine->GetNumberOfTimesteps() != 0))
	{
		throw std::logic_error("CartesianSubgrid: set native gather budget before registration/advance");
	}
	m_FieldGatherBudget = bytes;
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "operator_cartesiansubgrid.h"
#include "CSPrimBox.h"
#include "CSPropDumpBox.h"
#include "CSPropExcitation.h"
#include "CSPropProbeBox.h"
#include "engine.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/operator_ext_tfsf.h"
#include "tools/constants.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
unsigned int BaseLineIndex(const double *lines, unsigned int count, double value)
{
	if (!std::isfinite(value))
	{
		throw std::invalid_argument("CartesianSubgrid: finite refinement coordinates required");
	}
	const double *next = std::lower_bound(lines, lines + count, value);
	unsigned int index = next == lines + count ? count - 1 : unsigned(next - lines);
	if (index > 0 && std::abs(lines[index - 1] - value) < std::abs(lines[index] - value))
	{
		--index;
	}
	const double step = index == 0 ? lines[1] - lines[0]
						: index + 1 == count
							? lines[index] - lines[index - 1]
							: std::min(lines[index] - lines[index - 1], lines[index + 1] - lines[index]);
	const double roundoff =
		32 * std::numeric_limits<double>::epsilon() * std::max({std::abs(value), std::abs(lines[index]), step});
	// Permit only coordinate roundoff, additionally bounded by local spacing.
	// This is not geometric snapping to the nearest line.
	if (std::abs(value - lines[index]) > std::min(roundoff, step * 1e-8))
	{
		throw std::invalid_argument("CartesianSubgrid: box faces must coincide with exact base-grid lines");
	}
	return index;
}

//! Serial reference execution; field sampling does not allocate a full lattice.
class Engine_CartesianSubgrid : public Engine
{
  public:
	explicit Engine_CartesianSubgrid(const Operator_CartesianSubgrid *op) : Engine(op), m_Operator(op) {}
	bool IterateTS(unsigned int steps) override
	{
		Excitation *signal = m_Operator->Signal();
		for (unsigned int step = 0; step < steps; ++step)
		{
			const double drive = numTS < signal->GetLength() ? signal->GetVoltageSignal()[numTS] : 0;
			m_Operator->Grid().Step(m_Operator->GetTimestep(), drive);
			++numTS;
		}
		return true;
	}
	FDTD_FLOAT GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const override
	{
		const unsigned int pos[3] = {x, y, z};
		return GetVolt(n, pos);
	}
	FDTD_FLOAT GetVolt(unsigned int n, const unsigned int pos[3]) const override
	{
		if (pos[n] + 1 >= m_Operator->GetNumberOfLines(n))
		{
			return 0;
		}
		return m_Operator->Grid().Electric(n, m_Operator->Position(n, pos, false)) * m_Operator->GetEdgeLength(n, pos);
	}
	FDTD_FLOAT GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const override
	{
		const unsigned int pos[3] = {x, y, z};
		return GetCurr(n, pos);
	}
	FDTD_FLOAT GetCurr(unsigned int n, const unsigned int pos[3]) const override
	{
		if (pos[(n + 1) % 3] + 1 >= m_Operator->GetNumberOfLines((n + 1) % 3) ||
			pos[(n + 2) % 3] + 1 >= m_Operator->GetNumberOfLines((n + 2) % 3))
		{
			return 0;
		}
		return m_Operator->Grid().Magnetic(n, m_Operator->Position(n, pos, true)) *
			   m_Operator->GetEdgeLength(n, pos, true);
	}
	bool CalcFastEnergy(const unsigned int[3], double &electric, double &magnetic) const override
	{
		// Engine_Interface_FDTD multiplies by vacuum epsilon/mu. Return the
		// correctly integrated energy, not a count-dependent coarse sample sum.
		electric = 2 * m_Operator->Grid().ElectricEnergy() / EPS0;
		magnetic = 2 * m_Operator->Grid().MagneticEnergy() / MUE0;
		return true;
	}
	void SetVolt(unsigned int, unsigned int, unsigned int, unsigned int, FDTD_FLOAT) override
	{
		RejectWrite();
	}
	void SetVolt(unsigned int, const unsigned int[3], FDTD_FLOAT) override
	{
		RejectWrite();
	}
	void SetCurr(unsigned int, unsigned int, unsigned int, unsigned int, FDTD_FLOAT) override
	{
		RejectWrite();
	}
	void SetCurr(unsigned int, const unsigned int[3], FDTD_FLOAT) override
	{
		RejectWrite();
	}

  private:
	void RejectWrite() const
	{
		throw std::logic_error("CartesianSubgrid: use the core field/source API, not rectangular field writes");
	}
	const Operator_CartesianSubgrid *m_Operator;
};
} // namespace

Operator_CartesianSubgrid::Operator_CartesianSubgrid(const std::vector<std::vector<double>> &boxes,
													 const std::vector<CartesianSubgrid::Index> &ratios)
	: m_Boxes(boxes), m_Ratios(ratios)
{
}

Operator_CartesianSubgrid *Operator_CartesianSubgrid::New(const std::vector<double> &box, CartesianSubgrid::Index ratio)
{
	return NewRegions({box}, {ratio});
}

void Operator_CartesianSubgrid::ValidateRegions(const std::vector<std::vector<double>> &boxes,
												const std::vector<CartesianSubgrid::Index> &ratios)
{
	if (boxes.empty() || boxes.size() > 64 || boxes.size() != ratios.size())
	{
		throw std::invalid_argument("CartesianSubgrid: expected 1 through 64 boxes and matching axis ratios");
	}
	for (std::size_t i = 0; i < boxes.size(); ++i)
	{
		const auto &box = boxes[i];
		if (box.size() != 6 || !std::all_of(box.begin(), box.end(), [](double v) { return std::isfinite(v); }))
		{
			throw std::invalid_argument("CartesianSubgrid: expected six finite box coordinates: x0,y0,z0,x1,y1,z1");
		}
		CartesianSubgrid::ValidateRatio(ratios[i]);
		for (unsigned int n = 0; n < 3; ++n)
		{
			if (box[n] >= box[n + 3])
			{
				throw std::invalid_argument("CartesianSubgrid: box start must be below stop on every axis");
			}
		}
		for (std::size_t j = 0; j < i; ++j)
		{
			bool separated = false;
			for (unsigned int n = 0; n < 3; ++n)
			{
				separated = separated || box[n + 3] < boxes[j][n] || boxes[j][n + 3] < box[n];
			}
			if (!separated)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: regions must be separated; touching, overlap and nesting are unsupported");
			}
		}
	}
}

Operator_CartesianSubgrid *Operator_CartesianSubgrid::NewRegions(const std::vector<std::vector<double>> &boxes,
																 const std::vector<CartesianSubgrid::Index> &ratios)
{
	ValidateRegions(boxes, ratios);
	std::unique_ptr<Operator_CartesianSubgrid> op(new Operator_CartesianSubgrid(boxes, ratios));
	op->Init();
	return op.release();
}

void Operator_CartesianSubgrid::SetEngineType(Engine_Types::EngineType type)
{
	if (type != Engine::BASIC)
	{
		throw std::invalid_argument("CartesianSubgrid: the XYZ reference implementation requires engine=basic");
	}
	Operator::SetEngineType(type);
}

void Operator_CartesianSubgrid::SetBoundaryCondition(int *conditions)
{
	for (int i = 0; i < 6; ++i)
	{
		if (conditions[i] != 0)
		{
			throw std::invalid_argument("CartesianSubgrid: only outer PEC boundaries are supported");
		}
	}
	Operator::SetBoundaryCondition(conditions);
}

void Operator_CartesianSubgrid::AddExtension(Operator_Extension *extension)
{
	std::unique_ptr<Operator_Extension> owner(extension);
	if (!dynamic_cast<Operator_Ext_Excitation *>(extension) && !dynamic_cast<Operator_Ext_TFSF *>(extension))
	{
		throw std::invalid_argument("CartesianSubgrid: unsupported extension " + extension->GetExtensionName());
	}
	// Soft E sources are prepared on compact master edges below. Actual TFSF
	// primitives are rejected; SetupFDTD adds its extension unconditionally.
}

bool Operator_CartesianSubgrid::SetGeometryCSX(ContinuousStructure *geometry)
{
	if (!geometry || geometry->GetGrid()->GetMeshType() != CARTESIAN)
	{
		throw std::invalid_argument("CartesianSubgrid: Cartesian CSX geometry required");
	}
	if (!Operator::SetGeometryCSX(geometry))
	{
		return false;
	}
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i < numLines[n]; ++i)
		{
			if (!std::isfinite(discLines[n][i]) || (i > 0 && !(discLines[n][i] > discLines[n][i - 1])))
			{
				throw std::invalid_argument("CartesianSubgrid: strictly increasing finite base coordinates required");
			}
		}
	}
	ValidateProperties();
	return true;
}

void Operator_CartesianSubgrid::ValidateProperties() const
{
	if (GetBackgroundKappa() != 0 || GetBackgroundSigma() != 0)
	{
		throw std::invalid_argument("CartesianSubgrid: lossy background is not supported");
	}
	for (CSProperties *property : CSX->GetPropertyByType(CSProperties::ANY))
	{
		const int type = property->GetType();
		if (type & ~(CSProperties::EXCITATION | CSProperties::PROBEBOX | CSProperties::DUMPBOX))
		{
			throw std::invalid_argument(
				"CartesianSubgrid: only homogeneous background, soft E sources and E/H observations are supported");
		}
		if (auto *source = dynamic_cast<CSPropExcitation *>(property))
		{
			if (source->GetExcitType() != 0 || source->GetDelay() != 0)
			{
				throw std::invalid_argument("CartesianSubgrid: only zero-delay soft E sources are supported");
			}
		}
		if (auto *dump = dynamic_cast<CSPropDumpBox *>(property))
		{
			const int kind = dump->GetDumpType();
			if (kind != 0 && kind != 1 && kind != 10 && kind != 11)
			{
				throw std::invalid_argument("CartesianSubgrid: only E/H field dumps are supported");
			}
		}
		else if (auto *probe = dynamic_cast<CSPropProbeBox *>(property))
		{
			if (probe->GetProbeType() != 2 && probe->GetProbeType() != 3)
			{
				throw std::invalid_argument("CartesianSubgrid: only E/H probes are supported, not lumped ports");
			}
		}
		for (CSPrimitives *primitive : property->GetAllPrimitives())
		{
			if ((type & CSProperties::EXCITATION) && !dynamic_cast<CSPrimBox *>(primitive))
			{
				throw std::invalid_argument("CartesianSubgrid: soft E sources must use box primitives");
			}
			double box[6];
			if (primitive->HasTransform() || !primitive->GetBoundBox(box) ||
				!std::all_of(box, box + 6, [](double v) { return std::isfinite(v); }))
			{
				throw std::invalid_argument(
					"CartesianSubgrid: untransformed finite observation/source bounds required");
			}
		}
	}
}

int Operator_CartesianSubgrid::CalcECOperator(DebugFlags flags)
{
	if (flags != None)
	{
		throw std::invalid_argument("CartesianSubgrid: rectangular operator debug dumps are not supported");
	}
	CartesianSubgrid::Axes lines;
	for (unsigned int n = 0; n < 3; ++n)
	{
		for (unsigned int i = 0; i < numLines[n]; ++i)
		{
			lines[n].push_back((discLines[n][i] - discLines[n][0]) * gridDelta);
		}
	}
	std::vector<CartesianSubgrid::Region> regions;
	for (std::size_t i = 0; i < m_Boxes.size(); ++i)
	{
		CartesianSubgrid::Region region;
		region.ratio = m_Ratios[i];
		for (unsigned int n = 0; n < 3; ++n)
		{
			for (unsigned int side = 0; side < 2; ++side)
			{
				(side ? region.box.stop[n] : region.box.start[n]) =
					BaseLineIndex(discLines[n], numLines[n], m_Boxes[i][n + 3 * side]);
			}
		}
		regions.push_back(region);
	}
	m_Grid.reset(new CartesianSubgrid(lines, GetBackgroundEpsR() * EPS0, GetBackgroundMueR() * MUE0, regions));
	const double limit = m_Grid->StableTimestep();
	if (dT == 0)
	{
		dT = limit * m_TimeStepFactor;
	}
	if (!(dT > 0) || !std::isfinite(dT) || dT > limit)
	{
		throw std::invalid_argument("CartesianSubgrid: requested timestep exceeds the coupled bound");
	}
	if (m_Exc->GetSignalPeriod() != 0)
	{
		throw std::invalid_argument("CartesianSubgrid: periodic excitation is not yet supported");
	}
	m_Exc->Reset(dT);
	PrepareSources();
	return 0;
}

void Operator_CartesianSubgrid::PrepareSources()
{
	std::vector<std::pair<std::size_t, double>> sources;
	for (CSProperties *property : CSX->GetPropertyByType(CSProperties::EXCITATION))
	{
		auto *source = dynamic_cast<CSPropExcitation *>(property);
		if (!source || !source->GetEnabled())
		{
			continue;
		}
		std::size_t selected = 0;
		for (std::size_t edge = 0; edge < m_Grid->EdgeCount(); ++edge)
		{
			const unsigned int n = m_Grid->EdgeDirection(edge);
			if (!source->GetActiveDir(n) || source->GetExcitation(n) == 0)
			{
				continue;
			}
			double position[3];
			const auto physical = m_Grid->EdgePosition(edge);
			for (unsigned int k = 0; k < 3; ++k)
			{
				position[k] = physical[k] / gridDelta + discLines[k][0];
			}
			bool inside = false;
			for (CSPrimitives *primitive : source->GetAllPrimitives())
			{
				double box[6];
				primitive->GetBoundBox(box);
				bool hit = true;
				for (unsigned int k = 0; k < 3; ++k)
				{
					const double tolerance = 1e-8 * (discLines[k][1] - discLines[k][0]);
					hit = hit && position[k] >= box[2 * k] - tolerance && position[k] <= box[2 * k + 1] + tolerance;
				}
				inside = inside || hit;
				if (hit)
				{
					primitive->SetPrimitiveUsed(true);
				}
			}
			if (inside)
			{
				sources.emplace_back(edge, source->GetWeightedExcitation(n, position));
				++selected;
			}
		}
		if (selected == 0)
		{
			throw std::invalid_argument("CartesianSubgrid: source touches no unconstrained master edge");
		}
	}
	m_Grid->SetSources(sources);
}

Engine *Operator_CartesianSubgrid::CreateEngine()
{
	if (!m_Grid)
	{
		throw std::logic_error("CartesianSubgrid: build the operator first");
	}
	m_Engine = new Engine_CartesianSubgrid(this);
	return m_Engine;
}

CartesianSubgrid::Point Operator_CartesianSubgrid::Position(unsigned int n, const unsigned int pos[3], bool dual) const
{
	unsigned int copy[3] = {pos[0], pos[1], pos[2]};
	double coordinate[3];
	GetYeeCoords(n, copy, coordinate, dual);
	CartesianSubgrid::Point point;
	for (unsigned int k = 0; k < 3; ++k)
	{
		point[k] = std::clamp((coordinate[k] - discLines[k][0]) * gridDelta, 0.,
							  (discLines[k][numLines[k] - 1] - discLines[k][0]) * gridDelta);
	}
	return point;
}

double Operator_CartesianSubgrid::GetNumberCells() const
{
	return m_Grid ? m_Grid->CellCount() : 0;
}

void Operator_CartesianSubgrid::ShowStat() const
{
	std::cout << "Experimental compact Cartesian XYZ subgrid; " << m_Boxes.size()
			  << " region(s); serial host reference; " << m_Grid->CellCount() << " leaf cells; "
			  << m_Grid->StorageBytes() << " accounted bytes; full fine parent arrays: none" << std::endl;
	for (std::size_t i = 0; i < m_Ratios.size(); ++i)
	{
		std::cout << "  region " << i << ": ratio " << m_Ratios[i][0] << 'x' << m_Ratios[i][1] << 'x' << m_Ratios[i][2]
				  << std::endl;
	}
}

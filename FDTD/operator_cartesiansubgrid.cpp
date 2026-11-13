/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "operator_cartesiansubgrid.h"
#include "CSPrimBox.h"
#include "CSPrimCylinder.h"
#include "CSPrimCylindricalShell.h"
#include "CSPropDumpBox.h"
#include "CSPropExcitation.h"
#include "CSPropMaterial.h"
#include "CSPropLumpedElement.h"
#include "CSPropProbeBox.h"
#include "engine.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/operator_ext_tfsf.h"
#include "extensions/operator_ext_lumpedRLC.h"
#include "extensions/operator_ext_cartesiansubgrid_mur.h"
#include "tools/constants.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace
{
bool UnitWeight(const std::string &term)
{
	// CSXCAD stores a numeric weight without a term. Its value is checked
	// independently through the weighted material getters below.
	if (term.empty())
	{
		return true;
	}
	std::istringstream stream(term);
	double value = 0;
	if (!(stream >> value) || value != 1)
	{
		return false;
	}
	stream >> std::ws;
	return stream.eof();
}

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
	explicit Engine_CartesianSubgrid(const Operator_CartesianSubgrid *op) : Engine(op), m_Operator(op)
	{
		InitExtensions();
		SortExtensionByPriority();
	}
	bool IterateTS(unsigned int steps) override
	{
		Excitation *signal = m_Operator->Signal();
		for (unsigned int step = 0; step < steps; ++step)
		{
			const double drive = numTS < signal->GetLength() ? signal->GetVoltageSignal()[numTS] : 0;
			if (m_Eng_exts.empty())
			{
				m_Operator->Grid().Step(m_Operator->GetTimestep(), drive);
			}
			else
			{
				DoPreVoltageUpdates();
				m_Operator->Grid().UpdateElectric(m_Operator->GetTimestep(), drive);
				DoPostVoltageUpdates();
				Apply2Voltages();
				DoPreCurrentUpdates();
				m_Operator->Grid().UpdateMagnetic(m_Operator->GetTimestep());
				DoPostCurrentUpdates();
				Apply2Current();
			}
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
	int murAxis = -1;
	for (int i = 0; i < 6; ++i)
	{
		if (conditions[i] != 0 && conditions[i] != 1 && conditions[i] != 2)
		{
			throw std::invalid_argument("CartesianSubgrid: only outer PEC/PMC and scoped Mur are supported");
		}
		if (conditions[i] == 2)
		{
			if (murAxis != -1 && murAxis != i / 2)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: Mur faces must be parallel; adjacent Mur corners are unsupported");
			}
			murAxis = i / 2;
		}
	}
	Operator::SetBoundaryCondition(conditions);
}

void Operator_CartesianSubgrid::AddExtension(Operator_Extension *extension)
{
	std::unique_ptr<Operator_Extension> owner(extension);
	if (dynamic_cast<Operator_Ext_CartesianSubgrid_Mur *>(extension))
	{
		Operator::AddExtension(owner.release());
		return;
	}
	if (!dynamic_cast<Operator_Ext_Excitation *>(extension) &&
		!dynamic_cast<Operator_Ext_TFSF *>(extension) && !dynamic_cast<Operator_Ext_LumpedRLC *>(extension))
	{
		throw std::invalid_argument("CartesianSubgrid: unsupported extension " +
									extension->GetExtensionName());
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
		if (type != CSProperties::MATERIAL && type != CSProperties::METAL &&
			type != CSProperties::LUMPED_ELEMENT &&
			(type & ~(CSProperties::EXCITATION | CSProperties::PROBEBOX | CSProperties::DUMPBOX)))
		{
			throw std::invalid_argument(
				"CartesianSubgrid: only lossless isotropic materials, PEC, soft E sources "
				"and E/H observations are supported");
		}
		if (type == CSProperties::LUMPED_ELEMENT)
		{
			const auto *resistor = static_cast<CSPropLumpedElement *>(property);
			if (!(resistor->GetResistance() > 0) || !std::isfinite(resistor->GetResistance()) ||
				!std::isnan(resistor->GetInductance()) || !std::isnan(resistor->GetCapacity()) ||
				resistor->GetDirection() < 0 || resistor->GetDirection() > 2 ||
				resistor->GetLEtype() != CSPropLumpedElement::PARALLEL)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: lumped elements require a finite positive R only");
			}
		}
		if (type == CSProperties::MATERIAL)
		{
			auto *material = static_cast<CSPropMaterial *>(property);
			if (!material->GetIsotropy() || material->GetKappa() != 0 || material->GetSigma() != 0)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: material objects must be isotropic and lossless");
			}
			if (!(material->GetEpsilon() > 0) || !std::isfinite(material->GetEpsilon()) ||
				!(material->GetMue() > 0) || !std::isfinite(material->GetMue()))
			{
				throw std::invalid_argument(
					"CartesianSubgrid: positive finite material epsilon and permeability required");
			}
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (!UnitWeight(material->GetEpsilonWeightFunction(n)) ||
					!UnitWeight(material->GetMueWeightFunction(n)) ||
					!UnitWeight(material->GetKappaWeightFunction(n)) ||
					!UnitWeight(material->GetSigmaWeightFunction(n)))
				{
					throw std::invalid_argument(
						"CartesianSubgrid: spatial material weight functions are unsupported");
				}
				const double coordinate[3] = {1, 1, 1};
				if (material->GetEpsilonWeighted(n, coordinate) != material->GetEpsilon(n) ||
					material->GetMueWeighted(n, coordinate) != material->GetMue(n))
				{
					throw std::invalid_argument(
						"CartesianSubgrid: non-unit material weights are unsupported");
				}
			}
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
			if (probe->GetProbeType() < 0 || probe->GetProbeType() > 3)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: only voltage/current and E/H probes are supported");
			}
		}
		for (CSPrimitives *primitive : property->GetAllPrimitives())
		{
			if (type == CSProperties::LUMPED_ELEMENT && !dynamic_cast<CSPrimBox *>(primitive))
			{
				throw std::invalid_argument("CartesianSubgrid: resistors must use a line box");
			}
			if (type == CSProperties::METAL)
			{
				auto *cylinder = dynamic_cast<CSPrimCylinder *>(primitive);
				if (!dynamic_cast<CSPrimBox *>(primitive) && !cylinder)
				{
					throw std::invalid_argument(
						"CartesianSubgrid: PEC supports boxes, sheets and axis-aligned cylinders/shells");
				}
				if (cylinder)
				{
					unsigned int directions = 0;
					for (unsigned int n = 0; n < 3; ++n)
					{
						if (!std::isfinite(cylinder->GetCoord(2 * n)) ||
							!std::isfinite(cylinder->GetCoord(2 * n + 1)))
						{
							throw std::invalid_argument(
								"CartesianSubgrid: finite PEC cylinder coordinates required");
						}
						directions += cylinder->GetCoord(2 * n) != cylinder->GetCoord(2 * n + 1);
					}
					if (directions != 1 || !(cylinder->GetRadius() > 0) ||
						!std::isfinite(cylinder->GetRadius()))
					{
						throw std::invalid_argument("CartesianSubgrid: PEC cylinders require one finite "
													"nonzero axis and positive radius");
					}
					if (auto *shell = dynamic_cast<CSPrimCylindricalShell *>(cylinder))
					{
						if (!(shell->GetShellWidth() > 0) || !std::isfinite(shell->GetShellWidth()) ||
							shell->GetShellWidth() >= 2 * cylinder->GetRadius())
						{
							throw std::invalid_argument(
								"CartesianSubgrid: PEC shell requires positive wall and bore dimensions");
						}
					}
				}
			}
			if (type == CSProperties::MATERIAL && !dynamic_cast<CSPrimBox *>(primitive))
			{
				throw std::invalid_argument(
					"CartesianSubgrid: material objects must use axis-aligned volume boxes");
			}
			if ((type == CSProperties::MATERIAL || type == CSProperties::METAL) &&
				primitive->GetCoordinateSystem() != UNDEFINED_CS &&
				primitive->GetCoordinateSystem() != CARTESIAN)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: material/PEC primitives require Cartesian coordinates");
			}
			if ((type & CSProperties::EXCITATION) && !dynamic_cast<CSPrimBox *>(primitive))
			{
				throw std::invalid_argument("CartesianSubgrid: soft E sources must use box primitives");
			}
			double box[6];
			if (primitive->HasTransform() || !primitive->GetBoundBox(box) ||
				!std::all_of(box, box + 6, [](double v) { return std::isfinite(v); }))
			{
				throw std::invalid_argument(
					"CartesianSubgrid: untransformed finite primitive bounds required");
			}
			if (type == CSProperties::MATERIAL)
			{
				for (unsigned int n = 0; n < 3; ++n)
				{
					if (!(box[2 * n] < box[2 * n + 1]))
					{
						throw std::invalid_argument(
							"CartesianSubgrid: material boxes must have positive volume");
					}
				}
			}
			if (type == CSProperties::METAL && dynamic_cast<CSPrimBox *>(primitive))
			{
				unsigned int zero = 0;
				for (unsigned int n = 0; n < 3; ++n)
				{
					zero += box[2 * n] == box[2 * n + 1];
				}
				if (zero > 1)
				{
					throw std::invalid_argument(
						"CartesianSubgrid: PEC boxes require volume or one planar sheet, not wires/points");
				}
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
	// The primitive-list lookup stops at its first hit; CSXCAD requires this
	// list to be sorted by descending priority, not property insertion order.
	const auto primitives = CSX->GetAllPrimitives(true, CSProperties::MATERIAL);
	const auto geometry = CSX->GetAllPrimitives(
		true, (CSProperties::PropertyType)(CSProperties::MATERIAL | CSProperties::METAL));
	std::vector<std::array<double, 6>> bounds;
	std::vector<std::pair<unsigned int, std::array<double, 6>>> cylinders;
	for (CSPrimitives *primitive : geometry)
	{
		if (!dynamic_cast<CSPrimBox *>(primitive))
		{
			auto *cylinder = dynamic_cast<CSPrimCylinder *>(primitive);
			if (!cylinder)
			{
				throw std::invalid_argument("CartesianSubgrid: unsupported geometry added after validation");
			}
			std::array<double, 6> box;
			cylinder->GetBoundBox(box.data());
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (cylinder->GetCoord(2 * n) != cylinder->GetCoord(2 * n + 1))
				{
					cylinders.push_back({n, box});
				}
			}
			continue;
		}
		std::array<double, 6> box;
		primitive->GetBoundBox(box.data());
		bounds.push_back(box);
	}
	auto sample = [&](const CartesianSubgrid::Point &start, const CartesianSubgrid::Point &stop)
	{
		double lo[3], hi[3], centre[3], tolerance[3];
		for (unsigned int n = 0; n < 3; ++n)
		{
			lo[n] = start[n] / gridDelta + discLines[n][0];
			hi[n] = stop[n] / gridDelta + discLines[n][0];
			centre[n] = (lo[n] + hi[n]) / 2;
			if (!(hi[n] > lo[n]) || !(centre[n] > lo[n]) || !(centre[n] < hi[n]) ||
				!std::isfinite(centre[n]))
			{
				throw std::invalid_argument(
					"CartesianSubgrid: leaf material sampling is not representable in drawing coordinates");
			}
			tolerance[n] = std::min((hi[n] - lo[n]) * 1e-8,
									32 * std::numeric_limits<double>::epsilon() *
										std::max({std::abs(lo[n]), std::abs(hi[n]), hi[n] - lo[n]}));
		}
		// Centre sampling is exact only when every box is constant throughout
		// this leaf. Do not silently replace a subcell interface by its centre.
		for (const auto &box : bounds)
		{
			bool intersects = true, contains = true;
			int normal = -1;
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (box[2 * n] == box[2 * n + 1])
				{
					normal = int(n);
				}
			}
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (int(n) == normal)
				{
					intersects = intersects && box[2 * n] >= lo[n] - tolerance[n] &&
								 box[2 * n] <= hi[n] + tolerance[n];
					contains = contains && (std::abs(box[2 * n] - lo[n]) <= tolerance[n] ||
											std::abs(box[2 * n] - hi[n]) <= tolerance[n]);
					continue;
				}
				intersects =
					intersects && box[2 * n] < hi[n] - tolerance[n] && box[2 * n + 1] > lo[n] + tolerance[n];
				contains =
					contains && box[2 * n] <= lo[n] + tolerance[n] && box[2 * n + 1] >= hi[n] - tolerance[n];
			}
			if (intersects && !contains)
			{
				throw std::invalid_argument(
					"CartesianSubgrid: material/PEC box boundary cuts a leaf; align box "
					"faces to local cell boundaries");
			}
		}
		for (const auto &cylinder : cylinders)
		{
			const unsigned int axis = cylinder.first;
			const auto &box = cylinder.second;
			bool overlaps = true;
			for (unsigned int n = 0; n < 3; ++n)
			{
				if (n != axis)
				{
					overlaps = overlaps && box[2 * n] < hi[n] - tolerance[n] &&
							   box[2 * n + 1] > lo[n] + tolerance[n];
				}
			}
			for (unsigned int side = 0; side < 2; ++side)
			{
				const double end = box[2 * axis + side];
				if (overlaps && end > lo[axis] + tolerance[axis] && end < hi[axis] - tolerance[axis])
				{
					throw std::invalid_argument("CartesianSubgrid: PEC cylinder end cuts a leaf; align axial "
												"contacts to local boundaries");
				}
			}
		}
		auto *property = CSX->GetPropertyByCoordPriority(centre, primitives, true);
		auto *material = property ? static_cast<CSPropMaterial *>(property) : nullptr;
		return CartesianSubgrid::Medium{(material ? material->GetEpsilon() : GetBackgroundEpsR()) * EPS0,
										(material ? material->GetMue() : GetBackgroundMueR()) * MUE0};
	};
	auto pec = [&](unsigned int, const CartesianSubgrid::Point &start, const CartesianSubgrid::Point &stop)
	{
		double centre[3];
		for (unsigned int n = 0; n < 3; ++n)
		{
			centre[n] = (start[n] + stop[n]) / (2 * gridDelta) + discLines[n][0];
		}
		auto *property = CSX->GetPropertyByCoordPriority(centre, geometry, true);
		return property && property->GetType() == CSProperties::METAL;
	};
	const bool hasMetal = CSX->GetQtyPrimitives(CSProperties::METAL) != 0;
	CartesianSubgrid::BoundaryConditions boundaries;
	for (unsigned int i = 0; i < boundaries.size(); ++i)
	{
		boundaries[i] =
			m_BC[i] == 0 ? CartesianSubgrid::BoundaryType::PEC : CartesianSubgrid::BoundaryType::PMC;
	}
	m_Grid.reset(new CartesianSubgrid(
		lines, GetBackgroundEpsR() * EPS0, GetBackgroundMueR() * MUE0, regions, false, sample,
		hasMetal ? CartesianSubgrid::PECSelector(pec) : CartesianSubgrid::PECSelector(), boundaries));
	if (hasMetal && m_Grid->PECTraceCount() == 0)
	{
		throw std::invalid_argument(
			"CartesianSubgrid: no PEC edge resolved; check geometry, priorities and mesh");
	}
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
	PrepareResistors();
	PrepareSources();
	for (auto *extension : m_Op_exts)
	{
		if (!extension->BuildExtension())
		{
			throw std::invalid_argument("CartesianSubgrid: extension preparation failed");
		}
	}
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
			  << m_Grid->StorageBytes() << " accounted bytes; " << m_Grid->PECTraceCount()
			  << " PEC traces; full fine parent arrays: none" << std::endl;
	for (std::size_t i = 0; i < m_Ratios.size(); ++i)
	{
		std::cout << "  region " << i << ": ratio " << m_Ratios[i][0] << 'x' << m_Ratios[i][1] << 'x' << m_Ratios[i][2]
				  << std::endl;
	}
}

CartesianSubgrid::Point Operator_CartesianSubgrid::DrawingPoint(const double *drawing) const
{
	CartesianSubgrid::Point point;
	for (unsigned int n = 0; n < 3; ++n)
	{
		point[n] = (drawing[n] - discLines[n][0]) * gridDelta;
	}
	return point;
}

void Operator_CartesianSubgrid::PrepareResistors()
{
	std::vector<CartesianSubgrid::Resistor> ports;
	std::set<std::size_t> occupied;
	for (auto *property : CSX->GetPropertyByType(CSProperties::LUMPED_ELEMENT))
	{
		const auto *resistor = static_cast<CSPropLumpedElement *>(property);
		for (auto *primitive : property->GetAllPrimitives())
		{
			auto *box = static_cast<CSPrimBox *>(primitive);
			double start[3], stop[3];
			for (unsigned int n = 0; n < 3; ++n)
			{
				start[n] = box->GetCoord(2 * n);
				stop[n] = box->GetCoord(2 * n + 1);
				if (n != unsigned(resistor->GetDirection()) && start[n] != stop[n])
				{
					throw std::invalid_argument("CartesianSubgrid: native resistor currently requires one "
												"local edge, not a distributed box");
				}
			}
			const auto a = DrawingPoint(start), b = DrawingPoint(stop);
			if (!m_Grid->IsIndependentEdge(a, b))
			{
				throw std::invalid_argument("CartesianSubgrid: native resistor must span exactly one "
											"unconstrained independent local edge");
			}
			const auto path = m_Grid->VoltagePath(a, b);
			for (const auto &term : path.weights)
			{
				if (!occupied.insert(term.first).second)
				{
					throw std::invalid_argument(
						"CartesianSubgrid: overlapping native resistors are unsupported");
				}
			}
			ports.push_back({a, b, resistor->GetResistance(), 0});
			primitive->SetPrimitiveUsed(true);
		}
	}
	m_Grid->SetResistors(ports, dT);
}

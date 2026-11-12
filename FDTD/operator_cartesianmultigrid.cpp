/*
 * Based on operator_cylindermultigrid.cpp.
 * Copyright (C) 2010-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "operator_cartesianmultigrid.h"
#include "CSPropExcitation.h"
#include "cartesian_multigrid_coupling.h"
#include "engine.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/operator_ext_tfsf.h"
#include "tools/useful.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

Operator_CartesianMultiGrid::Operator_CartesianMultiGrid(double splitX)
	: Operator(), m_SplitX(splitX), m_Split(0), m_Coarse(NULL)
{
}

Operator_CartesianMultiGrid *Operator_CartesianMultiGrid::New(double splitX, unsigned int threads)
{
	if (!std::isfinite(splitX))
	{
		throw std::invalid_argument("CartesianMultiGrid: split X must be finite");
	}
	std::unique_ptr<Operator_CartesianMultiGrid> op(new Operator_CartesianMultiGrid(splitX));
	op->SetNumThreads(threads);
	op->Init();
	return op.release();
}

void Operator_CartesianMultiGrid::Init()
{
	Operator::Init();
	m_Coarse = Operator::New(m_numThreads);
}

Operator_CartesianMultiGrid::~Operator_CartesianMultiGrid()
{
	delete m_Coarse;
}

void Operator_CartesianMultiGrid::SetEngineType(Engine_Types::EngineType type)
{
	if (type != Engine::BASIC && type != Engine::SSE && type != Engine::CPU)
	{
		throw std::invalid_argument("CartesianMultiGrid: this prototype supports basic, SSE and CPU host engines only");
	}
	Operator::SetEngineType(type);
	m_Coarse->SetEngineType(type);
}

bool Operator_CartesianMultiGrid::SetupCSXGrid(CSRectGrid *grid)
{
	if (!Operator::SetupCSXGrid(grid))
	{
		return false;
	}
	if (numLines[1] < 9 || numLines[1] % 2 != 1)
	{
		throw std::invalid_argument("CartesianMultiGrid: Y must have an odd number of lines, at least nine");
	}
	// Initially qualify uniform tangential refinement only. Do not silently
	// interpret unrelated nonuniform intervals as an integer refinement ratio.
	const double dy = discLines[1][1] - discLines[1][0];
	for (unsigned int y = 1; y + 1 < numLines[1]; ++y)
	{
		if (std::abs(discLines[1][y + 1] - discLines[1][y] - dy) > 1e-9 * dy)
		{
			throw std::invalid_argument("CartesianMultiGrid: this prototype requires uniform Y spacing");
		}
	}
	m_Split = 0;
	for (unsigned int x = 0; x < numLines[0]; ++x)
	{
		if (std::abs(discLines[0][x] - m_SplitX) <= 1e-9 * GetMinMeshStep())
		{
			m_Split = x + 1;
			break;
		}
	}
	if (m_Split < 4 || m_Split > numLines[0] - 4)
	{
		throw std::invalid_argument(
			"CartesianMultiGrid: split must be an exact interior X line, with three collar cells per side");
	}
	return true;
}

bool Operator_CartesianMultiGrid::SetGeometryCSX(ContinuousStructure *geo)
{
	if (geo && geo->GetGrid()->GetMeshType() != CARTESIAN)
	{
		throw std::invalid_argument("CartesianMultiGrid: the CSX mesh must be Cartesian");
	}
	if (!Operator::SetGeometryCSX(geo))
	{
		return false;
	}
	CheckInterface();
	// Restore the shared CSX grid on every exit, including exceptions. Each
	// operator keeps its own copied lines; primitives remain in the same frame.
	struct RestoreGrid
	{
		CSRectGrid *grid;
		std::vector<double> lines[2];
		~RestoreGrid()
		{
			for (int axis = 0; axis < 2; ++axis)
			{
				grid->ClearLines(axis);
				for (double value : lines[axis])
				{
					grid->AddDiscLine(axis, value);
				}
			}
		}
	} restore;
	restore.grid = geo->GetGrid();
	for (int axis = 0; axis < 2; ++axis)
	{
		restore.lines[axis].assign(discLines[axis], discLines[axis] + numLines[axis]);
		restore.grid->ClearLines(axis);
	}
	for (unsigned int x = 0; x < m_Split; ++x)
	{
		restore.grid->AddDiscLine(0, discLines[0][x]);
	}
	for (unsigned int y = 0; y < numLines[1]; y += 2)
	{
		restore.grid->AddDiscLine(1, discLines[1][y]);
	}
	return m_Coarse->SetGeometryCSX(geo);
}

void Operator_CartesianMultiGrid::CheckInterface() const
{
	// This deliberately excludes material/conductor crossings, not just Yee
	// samples landing inside them. A later scheme must qualify such interfaces.
	const double lo = discLines[0][m_Split - 4];
	const double hi = discLines[0][m_Split + 2];
	for (CSProperties *property : CSX->GetPropertyByType(CSProperties::ANY))
	{
		const int supported = CSProperties::MATERIAL | CSProperties::METAL | CSProperties::EXCITATION |
							  CSProperties::PROBEBOX | CSProperties::DUMPBOX;
		if (property->GetType() & ~supported)
		{
			throw std::invalid_argument("CartesianMultiGrid: unsupported property " + property->GetName());
		}
		if (property->GetType() & CSProperties::EXCITATION)
		{
			CSPropExcitation *excitation = dynamic_cast<CSPropExcitation *>(property);
			if (excitation && excitation->GetExcitType() == 10)
			{
				throw std::invalid_argument("CartesianMultiGrid: TFSF sources are not supported");
			}
		}
		if (!(property->GetType() & (CSProperties::MATERIAL | CSProperties::METAL | CSProperties::EXCITATION |
									 CSProperties::PROBEBOX | CSProperties::DUMPBOX)))
		{
			continue;
		}
		for (size_t i = 0; i < property->GetQtyPrimitives(); ++i)
		{
			CSPrimitives *primitive = property->GetPrimitive(i);
			if (primitive->HasTransform())
			{
				throw std::invalid_argument("CartesianMultiGrid: transformed primitives are not qualified");
			}
			double bounds[6];
			const bool hasBounds = primitive->GetBoundBox(bounds);
			// CSPrimPoint reports UNDEFINED_CS for an implicit coordinate system;
			// its bounds are then in its Cartesian input mesh coordinates.
			const bool implicitCartesian = primitive->GetBoundBoxCoordSystem() == UNDEFINED_CS &&
										   primitive->GetCoordInputType() == CARTESIAN &&
										   primitive->GetCoordinateSystem() == UNDEFINED_CS;
			if (!hasBounds || (primitive->GetBoundBoxCoordSystem() != CARTESIAN && !implicitCartesian))
			{
				throw std::invalid_argument("CartesianMultiGrid: cannot verify a primitive's interface clearance");
			}
			for (double bound : bounds)
			{
				if (!std::isfinite(bound))
				{
					throw std::invalid_argument("CartesianMultiGrid: primitive bounds must be finite");
				}
			}
			if (bounds[0] <= hi && bounds[1] >= lo)
			{
				throw std::invalid_argument("CartesianMultiGrid: keep materials, conductors, sources and observations "
											"outside the interface collar");
			}
		}
	}
}

bool Operator_CartesianMultiGrid::GetYeeCoords(int component, unsigned int pos[3], double *coords, bool dual) const
{
	const bool inside = Operator::GetYeeCoords(component, pos, coords, dual);
	return inside && pos[0] >= m_Split - 1;
}

void Operator_CartesianMultiGrid::CalcStartStopLines(unsigned int &threads, std::vector<unsigned int> &start,
													 std::vector<unsigned int> &stop) const
{
	const auto jobs = AssignJobs2Threads(numLines[0] - m_Split + 1, threads, true);
	threads = jobs.size();
	start.resize(threads);
	stop.resize(threads);
	unsigned int begin = m_Split - 1;
	for (unsigned int i = 0; i < threads; ++i)
	{
		start[i] = begin;
		stop[i] = begin + jobs[i] - 1;
		begin = stop[i] + 1;
	}
}

void Operator_CartesianMultiGrid::SetBoundaryCondition(int *conditions)
{
	for (int face = 0; face < 6; ++face)
	{
		if (conditions[face] != 0)
		{
			throw std::invalid_argument(
				"CartesianMultiGrid: only PEC outer boundaries are qualified in this prototype");
		}
	}
	Operator::SetBoundaryCondition(conditions);
	m_Coarse->SetBoundaryCondition(conditions);
}

void Operator_CartesianMultiGrid::SetExcitationSignal(Excitation *excitation)
{
	Operator::SetExcitationSignal(excitation);
	m_Coarse->SetExcitationSignal(excitation);
}

void Operator_CartesianMultiGrid::SetMaterialStoreFlags(int type, bool value)
{
	Operator::SetMaterialStoreFlags(type, value);
	if (m_Coarse)
	{
		m_Coarse->SetMaterialStoreFlags(type, value);
	}
}

void Operator_CartesianMultiGrid::AddExtension(Operator_Extension *extension)
{
	std::unique_ptr<Operator_Extension> owner(extension);
	if (dynamic_cast<Operator_Ext_TFSF *>(extension))
	{
		// SetupFDTD adds this unconditionally. Actual TFSF primitives are rejected
		// by CheckInterface(), rather than silently losing their excitation.
		return;
	}
	if (!dynamic_cast<Operator_Ext_Excitation *>(extension))
	{
		throw std::invalid_argument("CartesianMultiGrid: unsupported extension " + extension->GetExtensionName());
	}
	std::unique_ptr<Operator_Extension> child(extension->Clone(m_Coarse));
	if (!child)
	{
		throw std::runtime_error("CartesianMultiGrid: failed to clone the excitation extension");
	}
	m_Coarse->AddExtension(child.release());
	Operator::AddExtension(owner.release());
}

int Operator_CartesianMultiGrid::CalcECOperator(DebugFlags flags)
{
	const double requested = dT;
	const double factor = m_TimeStepFactor;
	// First determine both native limits, then build both operators with one
	// identical step. No local subcycling, and no unverified forced larger step.
	dT = 0;
	m_TimeStepFactor = 1;
	m_Coarse->SetTimeStepMethod(m_TimeStepVar);
	m_Coarse->SetMaterialAvgMethod(m_MatAverageMethod);
	m_Coarse->SetTimestep(0);
	if (m_Coarse->CalcECOperator(flags) != 0 || Operator::CalcECOperator(flags) != 0)
	{
		return -1;
	}
	const double limit = std::min(dT, m_Coarse->GetTimestep());
	const double common = requested > 0 ? requested : limit * factor;
	if (!(common > 0) || !std::isfinite(common) || common > limit * (1 + 1e-12))
	{
		throw std::invalid_argument("CartesianMultiGrid: requested timestep exceeds a grid's stability bound");
	}
	dT = common;
	m_Coarse->SetTimestep(common);
	if (m_Coarse->CalcECOperator(flags) != 0 || Operator::CalcECOperator(flags) != 0)
	{
		return -1;
	}
	m_TimeStepFactor = factor;
	if (dT != m_Coarse->GetTimestep())
	{
		throw std::runtime_error("CartesianMultiGrid: the two operators disagree on the common timestep");
	}
	// The two native CFL limits alone do not bound a new coupled curl. For the
	// power-adjoint incidence use ||A||2 <= sqrt(||A||1 ||A||inf), where each
	// link of A has magnitude sqrt(VI_e*IV_h), already including the step.
	// Leapfrog requires ||A||2 < 2. Rebuild once if this conservative bound
	// tightens the common step. A user's explicitly forced unsafe step fails.
	m_CoupledStepNormBound = CoupledStepNormBound();
	if (m_CoupledStepNormBound >= 1.98)
	{
		if (requested > 0)
		{
			throw std::invalid_argument("CartesianMultiGrid: forced timestep exceeds the coupled-curl bound");
		}
		dT *= 1.98 / m_CoupledStepNormBound;
		m_TimeStepFactor = 1;
		m_Coarse->SetTimestep(dT);
		if (m_Coarse->CalcECOperator(flags) != 0 || Operator::CalcECOperator(flags) != 0)
		{
			return -1;
		}
		m_TimeStepFactor = factor;
		m_CoupledStepNormBound = CoupledStepNormBound();
		if (m_CoupledStepNormBound >= 2 || dT != m_Coarse->GetTimestep())
		{
			throw std::runtime_error("CartesianMultiGrid: coupled timestep qualification failed");
		}
	}
	SetupInterpolation();
	FillMissingMaterialStorage();
	return 0;
}

bool Operator_CartesianMultiGrid::Calc_ECPos(int component, const unsigned int *pos, double *circuit,
											 const std::vector<CSPrimitives *> &primitives) const
{
	// Ez is sampled by the coarse boundary restriction. Its interface mass is
	// the fine-side half dual cell plus, at even Y nodes, the coarse-side half
	// dual cell. Keeping a full fine dual cell at every node misplaces this
	// capacitance and creates an impedance discontinuity even in empty space.
	// Ey uses row-sum lumping of R^T Cc R, which preserves its constant mode.
	const bool result = Operator::Calc_ECPos(component, pos, circuit, primitives);
	if (result && component == 2 && pos[0] == m_Split - 1)
	{
		const unsigned int x = pos[0];
		const double left = (discLines[0][x] - discLines[0][x - 1]) / 2;
		const double right = (discLines[0][x + 1] - discLines[0][x]) / 2;
		const double massRatio = (right + (pos[1] % 2 == 0 ? 2 * left : 0)) / (left + right);
		circuit[0] *= massRatio;
		circuit[1] *= massRatio;
	}
	return result;
}

double Operator_CartesianMultiGrid::CoupledStepNormBound() const
{
	const Operator *grids[2] = {this, m_Coarse};
	std::vector<double> rows[2];
	unsigned int counts[2][3];
	for (int grid = 0; grid < 2; ++grid)
	{
		for (int axis = 0; axis < 3; ++axis)
		{
			counts[grid][axis] = grids[grid]->GetNumberOfLines(axis, true);
		}
		rows[grid].assign(size_t(3) * counts[grid][0] * counts[grid][1] * counts[grid][2], 0.);
	}
	double maxColumn = 0;
	for (int grid = 0; grid < 2; ++grid)
	{
		const Operator *op = grids[grid];
		for (unsigned int x = grid == 0 ? m_Split - 1 : 0; x + 1 < counts[grid][0]; ++x)
		{
			for (unsigned int y = 0; y + 1 < counts[grid][1]; ++y)
			{
				for (unsigned int z = 0; z + 1 < counts[grid][2]; ++z)
				{
					for (unsigned int component = 0; component < 3; ++component)
					{
						const double iv = op->GetIV(component, x, y, z);
						if (!(iv >= 0) || !std::isfinite(iv))
						{
							throw std::invalid_argument("CartesianMultiGrid: non-passive current coefficient");
						}
						double column = 0;
						auto link = [&](unsigned int electric, unsigned int ex, unsigned int ey, unsigned int ez)
						{
							const bool boundary = grid == 1 && ex == m_Split - 1 && electric != 0;
							const int target = boundary ? 0 : grid;
							const unsigned int firstY = boundary ? 2 * ey : ey;
							const unsigned int copies = boundary && electric == 1 ? 2 : 1;
							for (unsigned int offset = 0; offset < copies; ++offset)
							{
								const unsigned int yy = firstY + offset;
								const double vi = grids[target]->GetVI(electric, ex, yy, ez);
								if (!(vi >= 0) || !std::isfinite(vi))
								{
									throw std::invalid_argument("CartesianMultiGrid: non-passive voltage coefficient");
								}
								const double weight = std::sqrt(vi * iv);
								column += weight;
								const size_t index =
									((size_t(electric) * counts[target][0] + ex) * counts[target][1] + yy) *
										counts[target][2] +
									ez;
								rows[target][index] += weight;
							}
						};
						const unsigned int a = (component + 1) % 3;
						const unsigned int b = (component + 2) % 3;
						unsigned int pa[3] = {x, y, z};
						unsigned int pb[3] = {x, y, z};
						++pa[a];
						++pb[b];
						link(b, x, y, z);
						link(b, pa[0], pa[1], pa[2]);
						link(a, x, y, z);
						link(a, pb[0], pb[1], pb[2]);
						maxColumn = std::max(maxColumn, column);
					}
				}
			}
		}
	}
	double maxRow = 0;
	for (const auto &grid : rows)
	{
		maxRow = std::max(maxRow, *std::max_element(grid.begin(), grid.end()));
	}
	return std::sqrt(maxRow * maxColumn);
}

Engine *Operator_CartesianMultiGrid::CreateEngine()
{
	Operator::CreateEngine();
	Engine *coarse = NULL;
	try
	{
		coarse = m_Coarse->CreateEngine();
		m_Engine->SetSubGrid(coarse, new CartesianMultiGridCoupling(m_Engine, coarse, m_Interpolation, numLines));
	}
	catch (...)
	{
		if (m_Engine->GetSubGridEngine() != coarse)
		{
			delete coarse;
		}
		delete m_Engine;
		m_Engine = NULL;
		throw;
	}
	return m_Engine;
}

double Operator_CartesianMultiGrid::GetNumberCells() const
{
	return double(numLines[0] - m_Split + 1) * numLines[1] * numLines[2] + m_Coarse->GetNumberCells();
}

void Operator_CartesianMultiGrid::ShowStat() const
{
	std::cout << "Experimental Cartesian multi-grid: X < " << m_SplitX
			  << " uses every second Y line; common timestep; full parent storage retained" << std::endl;
	m_Coarse->ShowStat();
	Operator::ShowStat();
}

void Operator_CartesianMultiGrid::SetupInterpolation()
{
	// Processing reconstruction only. The actual interface uses voltage
	// restriction and its power-adjoint current transfer, not these weights.
	m_Interpolation = MultiGridInterpolation();
	m_Interpolation.split_pos = m_Split;
	for (int field = 0; field < 2; ++field)
	{
		for (int direction = 0; direction < 2; ++direction)
		{
			auto &first = field ? m_Interpolation.pos_i_2p[direction] : m_Interpolation.pos_v_2p[direction];
			auto &second = field ? m_Interpolation.pos_i_2pp[direction] : m_Interpolation.pos_v_2pp[direction];
			auto &wFirst = field ? m_Interpolation.w_i_2p[direction] : m_Interpolation.w_v_2p[direction];
			auto &wSecond = field ? m_Interpolation.w_i_2pp[direction] : m_Interpolation.w_v_2pp[direction];
			first.resize(numLines[1]);
			second.resize(numLines[1]);
			wFirst.assign(numLines[1], 0);
			wSecond.assign(numLines[1], 0);
			const bool dual = (field == 1) != (direction == 1);
			const unsigned int count = m_Coarse->GetNumberOfLines(1) - (dual ? 1 : 0);
			for (unsigned int y = 0; y < numLines[1]; ++y)
			{
				const double target = GetDiscLine(1, y, dual);
				unsigned int left = 0;
				while (left + 2 < count && m_Coarse->GetDiscLine(1, left + 1, dual) < target)
				{
					++left;
				}
				first[y] = left;
				second[y] = left + 1;
				if (field == 0 && direction == 1 && y + 1 == numLines[1])
				{
					continue; // no voltage edge beyond the positive PEC Y face
				}
				const double a = m_Coarse->GetDiscLine(1, left, dual);
				const double b = m_Coarse->GetDiscLine(1, left + 1, dual);
				double wa = (b - target) / (b - a);
				double wb = (target - a) / (b - a);
				if (direction == 1)
				{
					wa *= GetDiscDelta(1, y, field == 1) / m_Coarse->GetDiscDelta(1, left, field == 1);
					wb *= GetDiscDelta(1, y, field == 1) / m_Coarse->GetDiscDelta(1, left + 1, field == 1);
				}
				wFirst[y] = wa;
				wSecond[y] = wb;
			}
		}
	}
}

void Operator_CartesianMultiGrid::FillMissingMaterialStorage()
{
	unsigned int pos[3];
	double EffMat[4];
	for (int ny = 0; ny < 3; ++ny)
	{
		for (pos[0] = 0; pos[0] < m_Split - 1; ++pos[0])
		{
			for (pos[1] = 0; pos[1] < numLines[1]; ++pos[1])
			{
				std::vector<CSPrimitives *> vPrims =
					this->GetPrimitivesBoundBox(pos[0], pos[1], -1, CSProperties::MATERIAL);

				for (pos[2] = 0; pos[2] < numLines[2]; ++pos[2])
				{
					Calc_EffMatPos(ny, pos, EffMat, vPrims);

					if (m_epsR_ptr)
					{
						ArrayLib::ArrayNIJK<float> &m_epsR = *m_epsR_ptr;
						m_epsR(ny, pos[0], pos[1], pos[2]) = EffMat[0];
					}
					if (m_kappa_ptr)
					{
						ArrayLib::ArrayNIJK<float> &m_kappa = *m_kappa_ptr;
						m_kappa(ny, pos[0], pos[1], pos[2]) = EffMat[1];
					}
					if (m_mueR_ptr)
					{
						ArrayLib::ArrayNIJK<float> &m_mueR = *m_mueR_ptr;
						m_mueR(ny, pos[0], pos[1], pos[2]) = EffMat[2];
					}
					if (m_sigma_ptr)
					{
						ArrayLib::ArrayNIJK<float> &m_sigma = *m_sigma_ptr;
						m_sigma(ny, pos[0], pos[1], pos[2]) = EffMat[3];
					}
				}
			}
		}
	}
}

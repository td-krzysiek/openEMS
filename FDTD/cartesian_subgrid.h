/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_H
#define CARTESIAN_SUBGRID_H

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <utility>
#include <vector>

//! Experimental compact, one-level Cartesian XYZ refinement with a conforming edge trace.
/*! A refined base cell is replaced by ratio[0]*ratio[1]*ratio[2] leaves. Hanging tangential E
 *  values are constrained by the neighbouring coarse hexahedron's edge basis.
 *  The same prolongation is used in the curl and in the Galerkin electric
 *  mass, including its off-diagonal interface blocks. No cylindrical code,
 *  global fine field lattice, pointwise current interpolation or subcycling.
 *  Lossless isotropic media, constant within each leaf, and constrained PEC edges.
 */
class CartesianSubgrid
{
  public:
	using Point = std::array<double, 3>;
	using Index = std::array<unsigned int, 3>;
	using Axes = std::array<std::vector<double>, 3>;
	using Field = std::function<double(unsigned int, const Point &)>;
	enum class BoundaryType
	{
		PEC,
		PMC
	};
	//! X-/X+/Y-/Y+/Z-/Z+; value initialization retains PEC on every face.
	using BoundaryConditions = std::array<BoundaryType, 6>;
	struct Box
	{
		Index start;
		Index stop;
	};
	struct Region
	{
		Box box;
		Index ratio;
	};
	struct Leaf
	{
		std::array<int, 3> lo;
		std::array<int, 3> size;
	};
	struct Medium
	{
		double epsilon;
		double permeability;
	};
	//! Setup-only sampler: physical leaf bounds in metres, absolute SI coefficients.
	using MediumSampler = std::function<Medium(const Point &, const Point &)>;
	//! Setup-only PEC edge classifier in absolute SI coordinates.
	using PECSelector = std::function<bool(unsigned int, const Point &, const Point &)>;
	//! Setup-only exact line/contour weights on this grid's normalized state.
	struct Integral
	{
		bool magnetic = false;
		std::vector<std::pair<std::size_t, double>> weights;
	};
	struct Resistor
	{
		Point start, stop;
		double resistance;
		//! Norton current in amperes per unit Step sourceAmplitude (not a soft E kick).
		double current = 0;
	};
	Integral VoltagePath(const Point &start, const Point &stop) const;
	Integral CurrentContour(unsigned int normal, const Point &start, const Point &stop) const;
	double Evaluate(const Integral &integral) const;
	//! Exact independent local edge, used by the restricted native lumped adapter.
	bool IsIndependentEdge(const Point &start, const Point &stop) const;
	//! Dual contour bounds around a zero-width, node-centred native probe.
	std::pair<double, double> DualBounds(unsigned int axis, const Point &point) const;
	//! Fixed-step, passive midpoint resistor coupling; at most 64 small reference ports.
	void SetResistors(const std::vector<Resistor> &ports, double timestep);
	double ResistorVoltage(std::size_t port) const;
	double ResistorLoss() const
	{
		return m_ResistorLoss;
	}
	double SourceWork() const
	{
		return m_SourceWork;
	}

	CartesianSubgrid(Index cells, Point spacing, Box refined, double epsilon, double permeability,
					 bool uniformFine = false, Index ratio = {2, 2, 2});
	CartesianSubgrid(Index cells, Point spacing, double epsilon, double permeability,
					 const std::vector<Region> &regions, bool uniformFine = false);
	CartesianSubgrid(const Axes &lines, double epsilon, double permeability,
					 const std::vector<Region> &regions, bool uniformFine = false,
					 const MediumSampler &medium = MediumSampler(), const PECSelector &pec = PECSelector(),
					 const BoundaryConditions &boundaries = BoundaryConditions{});
	//! Actual local edge classification, including dependent interface traces.
	bool LeafEdgeIsPEC(std::size_t leaf, unsigned int axis, unsigned int first, unsigned int second) const;
	std::size_t PECTraceCount() const
	{
		return m_PECEdges.size();
	}
	Medium LeafMedium(std::size_t leaf) const
	{
		return m_Media.at(leaf);
	}
	//! Physical base-grid coordinates; integer Leaf extents index these intervals.
	const Axes &BaseLines() const
	{
		return m_Lines;
	}
	Point LeafStart(std::size_t leaf) const;
	Point LeafLengths(std::size_t leaf) const;
	static void ValidateRatio(const Index &ratio);
	//! Finest integer coordinate lattice; not the refinement of every region.
	Index LatticeRatio() const
	{
		return m_Ratio;
	}
	const std::vector<Region> &Regions() const
	{
		return m_Regions;
	}
	//! Compatibility alias for the single-region reference.

	Index Ratio() const
	{
		return m_Ratio;
	}
	double StableTimestep() const
	{
		return m_StableStep;
	}
	void Step(double timestep, double sourceAmplitude = 0);
	void SetSources(const std::vector<std::pair<std::size_t, double>> &sources);
	void SetElectric(const Field &field);
	void SetMagnetic(const Field &field);
	double Electric(unsigned int component, const Point &position) const;
	double Magnetic(unsigned int component, const Point &position) const;
	void AddElectric(std::size_t edge, double value);
	Point EdgePosition(std::size_t edge) const;
	unsigned int EdgeDirection(std::size_t edge) const;
	std::size_t EdgeCount() const
	{
		return m_Electric.size();
	}
	std::size_t FaceCount() const
	{
		return m_Magnetic.size();
	}
	std::size_t CellCount() const
	{
		return m_Leaves.size();
	}
	std::size_t StorageBytes() const;
	const std::vector<Leaf> &Leaves() const
	{
		return m_Leaves;
	}
	double Energy(double timestep = 0) const;
	double ElectricEnergy() const;
	double MagneticEnergy() const;
	//! Sum of the six oriented flux changes must vanish per leaf (D C = 0).
	double CurlDivergenceResidual(const Field &field) const;

  private:
	struct Edge
	{
		int axis;
		std::array<int, 3> lo;
		int size;
		bool operator<(const Edge &other) const;
	};
	struct Trace
	{
		std::vector<std::pair<std::size_t, double>> weights;
	};
	struct MassBlock
	{
		std::vector<std::size_t> edges;
		std::vector<double> lower;
		void Solve(std::vector<double> &values) const;
	};
	Edge LocalEdge(const Leaf &cell, unsigned int axis, unsigned int first, unsigned int second) const;
	bool Boundary(const Edge &edge) const;
	bool Hanging(const Edge &edge) const;
	const Box *InterfaceBox(const Edge &edge) const;
	Trace Restrict(const Edge &edge, bool ignoreLocalPEC = false) const;
	std::size_t RootIndex(const Index &cell) const;
	std::size_t FindLeaf(const Point &position) const;
	double TraceValue(const Trace &trace, const std::vector<double> &amplitude) const;
	void Build(const MediumSampler &medium, const PECSelector &pec);
	Point EdgePoint(const Edge &edge, double fraction) const;
	double FactorMass(const std::map<std::pair<std::size_t, std::size_t>, double> &offDiagonal);
	void Curl(const std::vector<double> &electric, std::vector<double> &magnetic) const;
	void CurlTranspose(const std::vector<double> &magnetic, std::vector<double> &electric) const;
	void SolveMass(std::vector<double> &values) const;
	void CheckStep(double timestep) const;
	Integral LineIntegral(bool magnetic, const Point &start, const Point &stop) const;
	void ApplyResistors(double timestep, double sourceAmplitude);
	double Coordinate(unsigned int axis, double lattice) const;
	double Length(const Leaf &leaf, unsigned int axis) const;

	Index m_Cells;
	Index m_Ratio;
	Axes m_Lines;
	std::vector<Region> m_Regions;
	BoundaryConditions m_Boundaries;
	double m_Epsilon;
	double m_Mu;
	bool m_UniformFine;
	double m_StableStep = 0;
	std::vector<Leaf> m_Leaves;
	std::vector<Medium> m_Media;
	std::vector<std::size_t> m_RootLeaves;
	std::vector<Index> m_RootRatios;
	std::map<Edge, std::size_t> m_EdgeMap;
	std::set<Edge> m_PECEdges;
	std::vector<unsigned char> m_PECFaces;
	std::vector<Edge> m_Edges;
	std::vector<double> m_SqrtMass, m_SqrtReluctance;
	std::vector<MassBlock> m_Blocks;
	std::vector<std::pair<std::pair<std::size_t, std::size_t>, double>> m_MassOffDiagonal;
	std::vector<std::size_t> m_RowStart, m_Columns;
	std::vector<double> m_Values;
	std::vector<double> m_Electric, m_Magnetic, m_Increment, m_Curl;
	std::vector<std::pair<std::size_t, double>> m_Sources;
	std::vector<Integral> m_PortVoltage, m_PortResponse;
	std::vector<double> m_PortConductance, m_PortCurrent, m_PortLower, m_PortRHS, m_PortMid, m_SourceCurl;
	double m_PortTimestep = 0, m_ResistorLoss = 0, m_SourceWork = 0;
};
#endif

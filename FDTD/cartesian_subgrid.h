/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_H
#define CARTESIAN_SUBGRID_H

#include "subgrid_curl.h"
#include "subgrid_edge_index.h"
#include "subgrid_leaf_index.h"
#include "gpu_cartesiansubgrid.h"
#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
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
    // Scoped absorbing extension needs the normalized curl and collar state;
    // the lossless core update itself is unchanged.
	friend class Operator_Ext_CartesianSubgrid_PML;
	friend class Engine_Ext_CartesianSubgrid_PML;
	friend class Operator_Ext_CartesianSubgrid_Conductivity;
	friend class Engine_Ext_CartesianSubgrid_Conductivity;
	// Domain preparation borrows normalized initial values only during its constructor.
	friend class CartesianSubgridDomainPlan;
  public:
	using Point = std::array<double, 3>;
	using Index = std::array<unsigned int, 3>;
	using Axes = std::array<std::vector<double>, 3>;
	using Field = std::function<double(unsigned int, const Point &)>;
	//! Internal storage policy: CPU reference patterns or spatial GPU-preparation blocks.
	enum class CurlLayout
	{
		Patterns,
		SpatialBlocks
	};
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
	//! Whole excluded conductor leaves: no bulk E/H energy or curl inside.
	using OpaqueSelector = std::function<bool(const Point &, const Point &)>;
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
	//! Connected axis-aligned whole local edges; no partial, interior-metal or PEC segments.
	bool IsLocalEdgePath(const Point &start, const Point &stop) const;
	//! Ordered oriented whole-edge segments; empty means an invalid geometric path.
	std::vector<std::pair<Point, Point>> LocalEdgePathSegments(const Point &start, const Point &stop) const;
	//! Geometry-normalized per-edge R for a resolved tensor-product native box.
	//! No hanging/missing edge, partial endpoint or opaque/PEC support is admitted.
	//! Returned independent conductances are consumed by the paired loss extension,
	//! not the bounded dense reference-port system in SetResistors.
	std::vector<Resistor> DistributedResistor(const Point &start, const Point &stop,
		unsigned int direction, double resistance) const;
	struct ResolvedResistor
	{
		Resistor resistor;
		std::size_t edge = 0, leaf = 0;
		unsigned int first = 0, second = 0;
	};
	//! Same geometry/normalization, retaining its already verified local support.
	//! Grid topology must remain unchanged while these setup descriptors are used.
	std::vector<ResolvedResistor> PrepareDistributedResistor(const Point &start, const Point &stop,
		unsigned int direction, double resistance) const;
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
					 const BoundaryConditions &boundaries = BoundaryConditions{},
					 const OpaqueSelector &opaque = OpaqueSelector(),
					 CurlLayout layout = CurlLayout::Patterns,
					 const Axes &childLines = Axes{}, std::size_t maxLeafCells = 1000000);
	//! Actual local edge classification, including dependent interface traces.
	bool LeafEdgeIsPEC(std::size_t leaf, unsigned int axis, unsigned int first, unsigned int second) const;
	std::size_t PECTraceCount() const
	{
		return m_PECEdges.size();
	}
	bool LeafIsOpaque(std::size_t leaf) const
	{
		return !m_OpaqueCells.empty() && m_OpaqueCells.at(leaf);
	}
	//! Direct local trace, avoiding global line lookup and quadrature scans.
	Integral LeafEdgeVoltage(std::size_t leaf, unsigned int axis, unsigned int first, unsigned int second) const;
	struct SurfaceFace
	{
		std::size_t leaf = 0;
		unsigned int normal = 0, side = 0;
	};
	//! Resolved exposed faces; unresolved coarse-face subdivision rejects.
	std::vector<SurfaceFace> OpaqueBoundaryFaces() const;
	Medium LeafMedium(std::size_t leaf) const
	{
		return m_Media.at(leaf);
	}
	//! Physical base-grid coordinates; integer Leaf extents index these intervals.
	const Axes &BaseLines() const
	{
		return m_Lines;
	}
	//! Setup/inspection lookup in physical coordinates, with the existing boundary tie policy.
	std::size_t LeafAt(const Point &point) const
	{
		return FindLeaf(point);
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
	//! Serial extension hooks may act between these ordered electric/magnetic sweeps.
	void UpdateElectric(double timestep, double sourceAmplitude = 0);
	void UpdateMagnetic(double timestep);
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
	std::size_t CurlStorageBytes() const
	{
		return m_CurlOperator.AllocationBytes();
	}
	std::size_t CurlPatternCount() const
	{
		return m_CurlOperator.PatternCount();
	}
	struct StructuredCurlInfo
	{
		std::size_t blocks = 0, regularRows = 0, exceptionRows = 0, bytes = 0;
	};
	StructuredCurlInfo StructuredCurlStorage() const
	{
		return {m_CurlOperator.RegularBlockCount(), m_CurlOperator.RegularRowCount(),
			m_CurlOperator.RowCount() - m_CurlOperator.RegularRowCount(), m_CurlOperator.AllocationBytes()};
	}
	//! C++ reference diagnostics; counts allocated factor/work entries, not process RSS.
	struct InterfaceStorageInfo
	{
		std::size_t blocks = 0, largestBlock = 0, largestHalfBandwidth = 0;
		std::size_t factorEntries = 0, denseEquivalentEntries = 0, workspaceEntries = 0;
		std::size_t inverseBoundSolves = 0, bipartiteBoundBlocks = 0;
	};
	InterfaceStorageInfo InterfaceStorage() const;
	//! Legacy diagnostic snapshot; materialized only when explicitly requested.
	//! Immutable shared storage keeps existing vector-reference callers valid.
	const std::vector<Leaf> &Leaves() const
	{
		auto cached = m_DiagnosticLeaves.Load();
		if (!cached)
		{
			std::shared_ptr<const std::vector<Leaf>> built =
				std::make_shared<std::vector<Leaf>>(m_Leaves.begin(), m_Leaves.end());
			if (m_DiagnosticLeaves.Publish(cached, built))
			{
				return *built;
			}
		}
		return *cached;
	}
	//! Compact diagnostics without materializing the legacy coordinate vector.
	const CartesianSubgridLeafIndex<Leaf> &LeafAddresses() const
	{
		return m_Leaves;
	}
	double Energy(double timestep = 0) const;
	double ElectricEnergy() const;
	double MagneticEnergy() const;
	//! Sum of the six oriented flux changes must vanish per leaf (D C = 0).
	double CurlDivergenceResidual(const Field &field) const;

	//! Checked upload and complete-step observation mirror for existing device backends.
	CartesianSubgridDeviceData DeviceData() const;
	//! Single-device preparation without the expanded curl/transpose oracle.
	CartesianSubgridDeviceData DeviceStateData() const;
	//! Checked domain metadata with global FTZ bounds, but no copied initial fields.
	CartesianSubgridDeviceData DeviceSetupMetadata() const;
	CartesianSubgridDeviceOperators DeviceCurls(bool spatialBlocks = true) const;
	//! Read one canonical normalized row into reusable setup-only storage.
	void CurlRow(std::size_t face, std::vector<std::pair<std::size_t, double>> &row) const;
	//! Exact ascending transpose row from local topology, without a global incoming CSR.
	void CurlTransposeRow(std::size_t edge, std::vector<std::pair<std::size_t, double>> &row) const;
	void ImportDeviceLedger(double resistorLoss, double sourceWork);
	void ImportDeviceFields(const float *electric, std::size_t electricCount,
							const float *magnetic, std::size_t magneticCount);

  private:
	CartesianSubgridDeviceData PrepareDeviceData(bool includeCurls, bool includeFields = true) const;
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
		std::size_t halfBandwidth = 0;
		//! Lower band, diagonal first in each row: L(i,j) at i*(b+1)+i-j.
		std::vector<double> lower;
		//! Solve a contiguous local RHS; caller gathers/scatters independent edges.
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
	void Build(const MediumSampler &medium, const PECSelector &pec, const OpaqueSelector &opaque, std::size_t expectedLeaves);
	void ConfigureRegularCurl();
	CartesianSubgridMetricView CurlMetrics() const
	{
		return {m_SqrtMass.data(), m_SqrtReluctance.data(), m_SqrtMass.size(), 6 * m_Leaves.size(), true};
	}
	Point EdgePoint(const Edge &edge, double fraction) const;
	double FactorMass(const std::map<std::pair<std::size_t, std::size_t>, double> &offDiagonal);
	void Curl(const std::vector<double> &electric, std::vector<double> &magnetic) const;
	void CurlTranspose(const std::vector<double> &magnetic, std::vector<double> &electric) const;
	void SolveMass(std::vector<double> &values);
	void CheckStep(double timestep) const;
	Integral LineIntegral(bool magnetic, const Point &start, const Point &stop) const;
	void ApplyResistors(double timestep, double sourceAmplitude);
	void ApplyResistors(double timestep, double sourceAmplitude,
		const std::vector<Integral> &responses, const std::vector<double> &lower);
	double Coordinate(unsigned int axis, double lattice) const;
	double LocalFraction(const Leaf &leaf, unsigned int axis, double position) const;
	double Length(const Leaf &leaf, unsigned int axis) const;

	Index m_Cells;
	Index m_Ratio;
	Axes m_Lines;
	//! Optional exact finest-lattice coordinates, including every base anchor.
	//! Empty axes retain the original equal-subdivision arithmetic.
	Axes m_ChildLines;
	std::vector<Region> m_Regions;
	BoundaryConditions m_Boundaries;
	double m_Epsilon;
	double m_Mu;
	bool m_UniformFine;
	CurlLayout m_CurlLayout;
	double m_StableStep = 0;
	CartesianSubgridLeafIndex<Leaf> m_Leaves;
	CartesianSubgridLeafSnapshot<Leaf> m_DiagnosticLeaves;
	std::vector<Medium> m_Media;
	std::vector<unsigned char> m_OpaqueCells;
	std::optional<CartesianSubgridEdgeIndex> m_EdgeIndex;
	std::set<Edge> m_PECEdges;
	std::vector<unsigned char> m_PECFaces;
	std::vector<Edge> m_Edges;
	//! Reluctance stores three metrics per leaf; six face states remain distinct.
	std::vector<double> m_SqrtMass, m_SqrtReluctance;
	std::vector<MassBlock> m_Blocks;
	std::vector<double> m_InterfaceWork;
	std::size_t m_InverseBoundSolves = 0, m_BipartiteBoundBlocks = 0;
	std::vector<std::pair<std::pair<std::size_t, std::size_t>, double>> m_MassOffDiagonal;
	CartesianSubgridCurl m_CurlOperator;
	std::vector<double> m_Electric, m_Magnetic, m_Increment, m_Curl;
	std::vector<std::pair<std::size_t, double>> m_Sources;
	std::vector<Integral> m_PortVoltage, m_PortResponse;
	std::vector<double> m_PortConductance, m_PortCurrent, m_PortLower, m_PortRHS, m_PortMid, m_SourceCurl;
	double m_PortTimestep = 0, m_ResistorLoss = 0, m_SourceWork = 0;
	double m_PendingMagneticTimestep = 0;
};
#endif

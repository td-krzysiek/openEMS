/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_H
#define CARTESIAN_SUBGRID_H

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <utility>
#include <vector>

//! Experimental compact, one-level Cartesian XYZ refinement with a conforming edge trace.
/*! A refined base cell is replaced by ratio[0]*ratio[1]*ratio[2] leaves. Hanging tangential E
 *  values are constrained by the neighbouring coarse hexahedron's edge basis.
 *  The same prolongation is used in the curl and in the Galerkin electric
 *  mass, including its off-diagonal interface blocks. No cylindrical code,
 *  global fine field lattice, pointwise current interpolation or subcycling.
 *  Homogeneous isotropic medium and outer PEC only in this numerical core.
 */
class CartesianSubgrid
{
  public:
	using Point = std::array<double, 3>;
	using Index = std::array<unsigned int, 3>;
	using Field = std::function<double(unsigned int, const Point &)>;
	struct Box
	{
		Index start;
		Index stop;
	};
	struct Leaf
	{
		std::array<int, 3> lo;
		std::array<int, 3> size;
	};

	CartesianSubgrid(Index cells, Point spacing, Box refined, double epsilon, double permeability,
					 bool uniformFine = false, Index ratio = {2, 2, 2});
	static void ValidateRatio(const Index &ratio);
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
	Trace Restrict(const Edge &edge) const;
	std::size_t RootIndex(const Index &cell) const;
	std::size_t FindLeaf(const Point &position) const;
	double TraceValue(const Trace &trace, const std::vector<double> &amplitude) const;
	void Build();
	double FactorMass(const std::map<std::pair<std::size_t, std::size_t>, double> &offDiagonal);
	void Curl(const std::vector<double> &electric, std::vector<double> &magnetic) const;
	void CurlTranspose(const std::vector<double> &magnetic, std::vector<double> &electric) const;
	void SolveMass(std::vector<double> &values) const;
	void CheckStep(double timestep) const;

	Index m_Cells;
	Index m_Ratio;
	Point m_Spacing;
	Box m_Box;
	double m_Epsilon;
	double m_Mu;
	bool m_UniformFine;
	double m_StableStep = 0;
	std::vector<Leaf> m_Leaves;
	std::vector<std::size_t> m_RootLeaves;
	std::map<Edge, std::size_t> m_EdgeMap;
	std::vector<Edge> m_Edges;
	std::vector<double> m_SqrtMass, m_SqrtReluctance;
	std::vector<MassBlock> m_Blocks;
	std::vector<std::pair<std::pair<std::size_t, std::size_t>, double>> m_MassOffDiagonal;
	std::vector<std::size_t> m_RowStart, m_Columns;
	std::vector<double> m_Values;
	std::vector<double> m_Electric, m_Magnetic, m_Increment, m_Curl;
	std::vector<std::pair<std::size_t, double>> m_Sources;
};
#endif

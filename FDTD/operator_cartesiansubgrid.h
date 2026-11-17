/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPERATOR_CARTESIANSUBGRID_H
#define OPERATOR_CARTESIANSUBGRID_H

#include "cartesian_subgrid.h"
#include "operator.h"
#include <memory>
#include <complex>
#include <optional>

class Operator_Ext_CartesianSubgrid_Conductivity;

class Engine_Field_Gather;
class Engine_Interface_FDTD;

//! Native API adapter for the experimental compact XYZ numerical core.
class Operator_CartesianSubgrid : public Operator
{
  public:
	static Operator_CartesianSubgrid *New(const std::vector<double> &box, CartesianSubgrid::Index ratio = {2, 2, 2});
	static Operator_CartesianSubgrid *NewRegions(const std::vector<std::vector<double>> &boxes,
												 const std::vector<CartesianSubgrid::Index> &ratios);
	static void ValidateRegions(const std::vector<std::vector<double>> &boxes,
								const std::vector<CartesianSubgrid::Index> &ratios);
	~Operator_CartesianSubgrid() override;
	void SetEngineType(Engine_Types::EngineType type) override;
	void SetChildLines(const CartesianSubgrid::Axes &drawingLines);
	//! Count-only opt-in budget; not a RAM estimate. Default remains one million.
	void SetMaxLeaves(std::size_t count);
	void SetMaterialSamples(unsigned int samples);
	//! Fixed execution placement, separate from physical geometry and XML.
	void SetDomainExecution(unsigned int axis, std::vector<double> drawingCuts,
		std::vector<std::unique_ptr<GPU_Backend>> backends);
	//! Startup registration; seal only after all native processing has been prepared.
	std::size_t RegisterObservation(const CartesianSubgrid::Integral &functional) const;
	void SealObservations();
	//! Startup-only shared gather budget; exhausted/late gathers retain the explicit full fallback.
	void SetFieldGatherBudget(std::size_t bytes);
	std::size_t FieldGatherBudget() const { return m_FieldGatherBudget; }
	Engine_Field_Gather *CreateFieldGather(const Engine_Interface_FDTD *interface, bool magnetic,
		const unsigned int count[3], unsigned int *const positions[3]) const;

	//! Native device DFT recipe IDs; -1 means an explicit unsupported/budget fallback.
	int CreateFieldDFT(bool magnetic, const std::vector<std::size_t> &observations, unsigned int frequencies) const;
	void AccumulateFieldDFT(int id, const std::vector<std::complex<float>> &weights) const;
	void ReadFieldDFT(int id, unsigned int frequency, std::vector<std::complex<double>> &sums) const;
	void SetFieldDFTBudget(std::size_t bytes);
	std::size_t FieldDFTStorageBytes() const { return m_FieldDFTUsed; }

	double EvaluateObservation(std::size_t id) const;
	struct ObservationInfo
	{
		std::size_t fullFieldReads = 0, selectedReadbacks = 0, energyReadbacks = 0;
		std::uint64_t energyDownloadBytes = 0;
		std::size_t hostObservationBytes = 0;
		bool sealed = false, device = false;
		std::size_t selectedElectricStates = 0, selectedMagneticStates = 0;
		std::uint64_t fieldDFTDownloadBytes = 0;

	};
	ObservationInfo NativeObservationInfo() const;
	//! Explicit complete mirror; refreshes lazily in sealed device execution.
	const CartesianSubgrid &ObservationGrid() const;
	void InvalidateObservations() const { m_ObservationsValid = false; }
	bool ObservationStateValid() const { return m_ObservationsValid; }
	//! Owner/processing thread only, after draining asynchronous readers and
	//! before publishing a paused witness. Returns logical worker acknowledgements.
	std::size_t FenceObservationEpoch(std::size_t expected) const;
	bool SetGeometryCSX(ContinuousStructure *geometry) override;
	void SetBoundaryCondition(int *conditions) override;
	void AddExtension(Operator_Extension *extension) override;
	int CalcECOperator(DebugFlags flags = None) override;
	Engine *CreateEngine() override;
	double GetNumberCells() const override;
	void ShowStat() const override;
	//! Raw setup/reference grid. After SealObservations use ObservationGrid() for live device fields.
	CartesianSubgrid &Grid() const
	{
		return *m_Grid;
	}
	CartesianSubgrid::Point Position(unsigned int n, const unsigned int pos[3], bool dual) const;
	CartesianSubgrid::Point DrawingPoint(const double *drawing) const;
	Excitation *Signal() const
	{
		return m_Exc;
	}

  private:
	Operator_CartesianSubgrid(const std::vector<std::vector<double>> &boxes,
							  const std::vector<CartesianSubgrid::Index> &ratios);
	void ValidateProperties() const;
	void PrepareSources();
	std::vector<CartesianSubgrid::ResolvedResistor> PrepareResistors();
	std::vector<std::vector<double>> m_Boxes;
	std::vector<CartesianSubgrid::Index> m_Ratios;
	CartesianSubgrid::Axes m_ChildDrawingLines;
	std::size_t m_MaxLeaves = 1000000;
	unsigned int m_MaterialSamples = 0;
	std::unique_ptr<CartesianSubgrid> m_Grid;
	unsigned int m_DomainAxis = 0;
	std::vector<double> m_DomainCuts;
	std::vector<std::unique_ptr<GPU_Backend>> m_DomainBackends;
	mutable bool m_ObservationsValid = true;
	mutable std::vector<CartesianSubgrid::Integral> m_ObservationFunctions;
	mutable std::vector<unsigned char> m_ObservationDeviceOnly;
	struct FieldDFTPlan
	{
		std::size_t deviceID = 0;
		unsigned int frequencies = 0;
		std::optional<unsigned int> lastEpoch;
		std::vector<CartesianSubgrid::Integral> samples;
	};
	mutable std::vector<FieldDFTPlan> m_FieldDFTPlans;
	std::size_t m_FieldDFTBudget = std::size_t(256) << 20;
	mutable std::size_t m_FieldDFTUsed = 0;

	bool m_ObservationPlanSealed = false;
	mutable std::size_t m_FieldGatherBudget = std::size_t(256) << 20, m_FieldGatherUsed = 0;
	//! Non-owning handle; the operator extension list retains ownership.
	Operator_Ext_CartesianSubgrid_Conductivity *m_NativeConductivity = nullptr;
};
#endif

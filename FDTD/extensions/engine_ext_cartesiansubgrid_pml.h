/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef ENGINE_EXT_CARTESIANSUBGRID_PML_H
#define ENGINE_EXT_CARTESIANSUBGRID_PML_H

#include "engine_extension.h"
#include "operator_ext_cartesiansubgrid_pml.h"
#include "FDTD/gpu_cartesiansubgrid_pml.h"

class Engine_Ext_CartesianSubgrid_Conductivity;

class Engine_Ext_CartesianSubgrid_PML : public Engine_Extension
{
  public:
	explicit Engine_Ext_CartesianSubgrid_PML(Operator_Ext_CartesianSubgrid_PML *op);
	void DoPreVoltageUpdates() override;
	void DoPostVoltageUpdates() override;
	void Apply2Voltages() override;
	void DoPreCurrentUpdates() override;
	void DoPostCurrentUpdates() override;
	void Apply2Current() override;
	//! Signed field work, not a nonnegative PML energy or radiated power.
	double BoundaryWorkOut() const { return m_WorkOut; }
	//! Fresh paired device state from already checked normalized operator rows.
	std::unique_ptr<GPU_CartesianSubgridPML> CreateDevice(
		GPU_Backend &backend, const GPU_CartesianSubgridFields &fields);
	void ReadDeviceDiagnostics(GPU_CartesianSubgridPML &device);
	void BindDomainDevice(const Engine_Ext_CartesianSubgrid_Conductivity *advance);
	void ReadDomainDiagnostics(GPU_CartesianSubgrid &device);
	//! Cold tables only; no running convolution history is migrated.
	CartesianSubgridPMLData PrepareDomainData(
		const Engine_Ext_CartesianSubgrid_Conductivity *advance = nullptr) const;
	const CartesianSubgrid &Grid() const;
	const Operator_Ext_CartesianSubgrid_PML &PreparedOperator() const { return *m_PML; }

  private:
	CartesianSubgridPMLData PrepareDeviceData() const;
	struct State
	{
		std::vector<double> history, correction, previous, target;
		explicit State(std::size_t size);
	};
	void Prepare(const std::vector<Operator_Ext_CartesianSubgrid_PML::Derivative> &nodes, State &state,
				 const std::vector<double> &source, const std::vector<double> &field, double sign);
	void Finish(const std::vector<Operator_Ext_CartesianSubgrid_PML::Derivative> &nodes, State &state,
				const std::vector<double> &field);
	void Apply(const std::vector<Operator_Ext_CartesianSubgrid_PML::Derivative> &nodes, const State &state,
			   std::vector<double> &field);
	Operator_Ext_CartesianSubgrid_PML *m_PML;
	State m_Electric, m_Magnetic;
	double m_WorkOut = 0;
	bool m_DeviceConfigured = false, m_HostAdvanced = false, m_DomainConfigured = false;
};
#endif

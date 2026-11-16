/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef ENGINE_EXT_CARTESIANSUBGRID_CONDUCTIVITY_H
#define ENGINE_EXT_CARTESIANSUBGRID_CONDUCTIVITY_H

#include "operator_ext_cartesiansubgrid_conductivity.h"
#include "FDTD/gpu_cartesiansubgrid_conductivity.h"

class Engine_Ext_CartesianSubgrid_PML;

class Engine_Ext_CartesianSubgrid_Conductivity : public Engine_Ext_CartesianSubgrid_Advance
{
  public:
	explicit Engine_Ext_CartesianSubgrid_Conductivity(Operator_Ext_CartesianSubgrid_Conductivity *op);
	void UpdateElectric(double timestep, double drive) override;
	void UpdateMagnetic(double timestep) override;
	std::unique_ptr<GPU_CartesianSubgrid> CreateDevice(GPU_Backend &backend) override;
	void ReadDeviceDiagnostics(GPU_CartesianSubgrid &device) override;
	double ElectricLoss() const { return m_ElectricLoss[0] + m_ElectricLoss[1]; }
	double MagneticLoss() const { return m_MagneticLoss[0] + m_MagneticLoss[1]; }
	double SurfaceLoss() const { return m_SurfaceLoss[0] + m_SurfaceLoss[1]; }
	double AdditionalStoredEnergy() const override;
	//! Cold-state C++ domain preparation; no live history migration is implied.
	CartesianSubgridConductivityPrepared PrepareDomainData(
		const Engine_Ext_CartesianSubgrid_PML *pml = nullptr) const;
	const Operator_Ext_CartesianSubgrid_Conductivity &PreparedOperator() const { return *m_Conductivity; }
	const CartesianSubgrid &Grid() const { return m_Grid; }

  private:
	static void Add(std::array<double, 2> &sum, double value);
	void Check(double timestep) const;
	CartesianSubgridConductivityPrepared PrepareDeviceData(bool expanded) const;
	Operator_Ext_CartesianSubgrid_Conductivity *m_Conductivity;
	CartesianSubgrid &m_Grid;
	std::vector<double> m_OldElectric, m_SourceCurl, m_Scratch;
	std::array<double, 2> m_ElectricLoss{}, m_MagneticLoss{};
	std::array<double, 2> m_SurfaceLoss{};
	std::vector<double> m_SurfaceState;
	double m_DeviceSurfaceEnergy = 0;
	bool m_UsingDevice = false;
	bool m_Advanced = false;
};
#endif

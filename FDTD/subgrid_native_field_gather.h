/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_NATIVE_FIELD_GATHER_H
#define SUBGRID_NATIVE_FIELD_GATHER_H

#include "engine_field_gather.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <vector>

class Operator_CartesianSubgrid;
class Engine;

//! Startup-only prepared SI samples, delivered through the existing raw-region interpolator.
//! The interface/operator/engine must outlive this gather; placement and modes remain immutable.
class Field_Gather_Subgrid : public Field_Gather_FDTD
{
  public:
	class BudgetExceeded : public std::runtime_error
	{
	  public:
		BudgetExceeded() : std::runtime_error("CartesianSubgrid: native gather table budget exceeded") {}
	};
	Field_Gather_Subgrid(const Engine_Interface_FDTD *interface, const Operator_CartesianSubgrid *op,
		bool magnetic, const unsigned int count[3], unsigned int *const positions[3], std::size_t budget);
	void Evaluate(std::size_t first, std::size_t last, ArrayLib::ArrayNIJK<float> &field,
		const float *source = nullptr) const override;
	//! Conservative table/registry estimate, not process RSS or total device storage.
	std::size_t EstimatedStorageBytes() const { return m_Estimate; }
	std::uint64_t Readbacks() const;

	int CreateDeviceDFT(const Engine_Interface_FDTD *caller, unsigned int frequencies) const;
	void AccumulateDeviceDFT(const Engine_Interface_FDTD *caller, const std::vector<std::complex<float>> &weights) const;
	void ReadDeviceDFT(const Engine_Interface_FDTD *caller, unsigned int frequency,
		ArrayLib::ArrayNIJK<std::complex<float>> &field) const;
	bool DFTOutputMatches(const Engine_Interface_FDTD *caller,
		const ArrayLib::ArrayNIJK<std::complex<float>> &field) const;
	unsigned int DFTFrequencyCount() const { return m_DFTFrequencies; }


  private:
	struct Sample
	{
		std::size_t observation = std::numeric_limits<std::size_t>::max();
		double length = 0;
	};
	void CheckInterface(const Engine_Interface_FDTD *caller) const;
	const Operator_CartesianSubgrid *m_Local;
	const Engine *m_EngineIdentity;
	mutable int m_DFTId = -1;
	mutable unsigned int m_DFTFrequencies = 0;
	std::vector<Sample> m_Samples;
	mutable std::vector<float> m_Raw;
	mutable std::mutex m_Mutex;
	mutable std::uint64_t m_Readbacks = 0;
	mutable unsigned int m_Epoch = std::numeric_limits<unsigned int>::max();
	std::size_t m_Estimate = 0, m_Lines = 0;
	int m_Mode = 0;
};
#endif

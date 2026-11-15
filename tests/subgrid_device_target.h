/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef TEST_SUBGRID_DEVICE_TARGET_H
#define TEST_SUBGRID_DEVICE_TARGET_H

#include "FDTD/gpu_backend.h"
#include <stdexcept>
#include <string>

//! A physical qualification request may not silently use another platform or CPU.
struct SubgridDeviceTarget
{
	std::string backend = "reference", expectedPrefix = "reference (";

	static SubgridDeviceTarget FromFlag(const std::string &flag)
	{
		if (flag == "--metal")
		{
			return {"auto", "Metal ("};
		}
		if (flag == "--cuda")
		{
#ifdef OPENEMS_SUBGRID_TEST_CUDA
			g_settings.SetGPUDevice(0);
			return {"auto", "HIP ("};
#else
			throw std::invalid_argument("CUDA qualification requested from a test built without NVIDIA/CUDA support");
#endif
		}
		throw std::invalid_argument("expected --metal or --cuda");
	}

	void Verify(const GPU_Backend &selected) const
	{
		if (selected.GetName().find(expectedPrefix) != 0)
		{
			throw std::runtime_error("subgrid qualification selected the wrong backend: " + selected.GetName());
		}
	}
};
#endif

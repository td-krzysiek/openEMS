/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2026 openEMS contributors
 */
#include "openems.h"
#include "Common/processing.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "tinyxml.h"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <sys/stat.h>

namespace
{
bool Exists(const char* path)
{
	struct stat info;
	return stat(path, &info) == 0;
}

bool ValidGateID(const std::string& id)
{
	if (id.empty() || id.size() > 64)
	{
		return false;
	}
	for (char c : id)
	{
		if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
			(c >= '0' && c <= '9') || c == '-' || c == '_'))
		{
			return false;
		}
	}
	return true;
}


}

void openEMS::ValidateStartupGate() const
{
	if (m_startupGateStep == 0)
	{
		if (!m_startupGateID.empty())
		{
			throw std::runtime_error("Startup gate: nonzero --startup-gate is required");
		}
		return;
	}
	if (!ValidGateID(m_startupGateID) || m_startupGateStep >= NrTS ||
		m_startupGateTimeout == 0 || m_startupGateTimeout > 3600)
	{
		throw std::runtime_error("Startup gate: invalid ID, timestep or timeout");
	}
	for (const char* path : {"startup-gate.xml", "startup-gate.xml.tmp",
		"STARTUP_CONTINUE", "STARTUP_REJECT"})
	{
		if (Exists(path))
		{
			throw std::runtime_error("Startup gate: stale control file; use a fresh directory");
		}
	}
}

bool openEMS::WaitStartupGate(double energy, double maxEnergy, double iterationSeconds)
{
	return WaitRunGate(m_startupGateID, "startup", 0, energy, maxEnergy, iterationSeconds);
}

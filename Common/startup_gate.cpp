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

bool HasDecision(const char* path, const std::string& id)
{
	if (!Exists(path))
	{
		return false;
	}
	std::ifstream file(path, std::ios::binary);
	std::string token;
	char c;
	while (token.size() <= 65 && file.get(c))
	{
		token += c;
	}
	if (!token.empty() && token.back() == '\n')
	{
		token.pop_back();
	}
	if (!token.empty() && token.back() == '\r')
	{
		token.pop_back();
	}
	if (token != id)
	{
		throw std::runtime_error("Startup gate: decision ID mismatch");
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
	// Flush existing samples only; do not call Process() or PostProcess() again.
	for (size_t n = 0; n < PA->GetNumberOfProcessings(); ++n)
	{
		Processing* processing = PA->GetProcessing(n);
		processing->FinishAsync();
		processing->FlushData();
	}
	TiXmlDocument document;
	TiXmlElement* header = new TiXmlElement("StartupGate");
	header->SetAttribute("id", m_startupGateID);
	header->SetAttribute("timestep", std::to_string(FDTD_Eng->GetNumberOfTimesteps()));
	header->SetDoubleAttribute("timestepSeconds", FDTD_Op->GetTimestep());
	header->SetDoubleAttribute("iterationSeconds", iterationSeconds);
	header->SetAttribute("maxTimesteps", std::to_string(NrTS));
	header->SetDoubleAttribute("endCriteria", endCrit);
	header->SetAttribute("converged", "false");
	document.LinkEndChild(header);
	if (!document.SaveFile("startup-gate.xml.tmp") ||
		std::rename("startup-gate.xml.tmp", "startup-gate.xml") != 0)
	{
		throw std::runtime_error("Startup gate: cannot publish witness");
	}
	std::cout << "Startup gate: waiting at timestep " << FDTD_Eng->GetNumberOfTimesteps()
		<< "; id=" << m_startupGateID << std::endl;
	const auto started = std::chrono::steady_clock::now();
	while (true)
	{
		if (CheckAbortCond())
		{
			return false;
		}
		if (m_dftSnapshots)
		{
			CheckDFTSnapshot(energy, maxEnergy);
		}
		if (HasDecision("STARTUP_REJECT", m_startupGateID))
		{
			throw std::runtime_error("Startup gate: qualification rejected");
		}
		if (HasDecision("STARTUP_CONTINUE", m_startupGateID))
		{
			std::cout << "Startup gate: accepted; continuing without reinitialization" << std::endl;
			return true;
		}
		if (std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
			>= m_startupGateTimeout)
		{
			throw std::runtime_error("Startup gate: qualification timed out");
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
}

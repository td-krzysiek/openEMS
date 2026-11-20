/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2026 openEMS contributors
 */
#include "openems.h"
#include "Common/processing.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "FDTD/excitation.h"
#include "tinyxml.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <sys/stat.h>

namespace
{
bool Exists(const std::string& path)
{
	struct stat info;
	return stat(path.c_str(), &info) == 0;
}

std::string Precise(double value)
{
	std::ostringstream text;
	text << std::setprecision(17) << value;
	return text.str();
}

bool ValidID(const std::string& value)
{
	if (value.empty() || value.size() > 48)
	{
		return false;
	}
	for (char c : value)
	{
		if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
			(c >= '0' && c <= '9') || c == '-' || c == '_'))
		{
			return false;
		}
	}
	return true;
}

bool Decision(const char* path, const std::string& id)
{
	std::ifstream file(path, std::ios::binary);
	if (!file)
	{
		return false;
	}
	std::string token;
	std::getline(file, token);
	if (token != id)
	{
		throw std::runtime_error(std::string(path).find("STARTUP_") == 0
			? "Startup gate: decision ID mismatch" : "Runtime gate: decision ID mismatch");
	}
	if (std::remove(path) != 0)
	{
		throw std::runtime_error("Runtime gate: cannot consume decision");
	}
	return true;
}
}

void openEMS::ValidateRuntimeGates()
{
	m_runtimeGateLevels.clear();
	m_runtimeGateIndex = 0;
	if (m_runtimeGateSchedule.empty() && !m_runtimeControl)
	{
		if (!m_runtimeGateID.empty())
		{
			throw std::runtime_error("Runtime gate: --decay-gates is required");
		}
		return;
	}
	if ((!m_runtimeGateSchedule.empty() && (!m_dftSnapshots || !m_exactEndCriteria)) || !ValidID(m_runtimeGateID) ||
		m_runtimeGateTimeout == 0 || m_runtimeGateTimeout > 3600)
	{
		throw std::runtime_error("Runtime gate: require snapshots, exact-endcriteria, valid ID and timeout 1..3600");
	}
	std::istringstream input(m_runtimeGateSchedule);
	std::string part;
	int previous = 0;
	while (std::getline(input, part, ','))
	{
		size_t length = 0;
		double value = std::stod(part, &length);
		if (length != part.size() || !std::isfinite(value) || value != std::floor(value) ||
			value >= previous || value < -300 || m_runtimeGateLevels.size() >= 32)
		{
			throw std::runtime_error("Runtime gate: thresholds must be decreasing integers in -1..-300 dB");
		}
		previous = static_cast<int>(value);
		m_runtimeGateLevels.push_back(previous);
	}
	if (!m_runtimeGateSchedule.empty() && (m_runtimeGateLevels.empty() || m_runtimeGateSchedule.back() == ','))
	{
		throw std::runtime_error("Runtime gate: invalid threshold list");
	}
	for (const char* path : {"runtime-gate.xml", "RUNTIME_PAUSE", "RUNTIME_CONTINUE", "RUNTIME_STOP"})
	{
		if (Exists(path))
		{
			throw std::runtime_error("Runtime gate: stale control file; use a fresh directory");
		}
	}
}

void openEMS::FlushProcessingForGate()
{
	for (size_t n = 0; n < PA->GetNumberOfProcessings(); ++n)
	{
		Processing* processing = PA->GetProcessing(n);
		processing->FinishAsync();
		processing->FlushData();
	}
}

bool openEMS::WaitRunGate(const std::string& id, const std::string& reason, int threshold,
	double energy, double maxEnergy, double iterationSeconds)
{
	const bool startup = reason == "startup";
	const unsigned int timeout = startup ? m_startupGateTimeout : m_runtimeGateTimeout;
	const char* witness = startup ? "startup-gate.xml" : "runtime-gate.xml";
	const char* temporary = startup ? "startup-gate.xml.tmp" : "runtime-gate.xml.tmp";
	const char* proceed = startup ? "STARTUP_CONTINUE" : "RUNTIME_CONTINUE";
	const char* stop = startup ? "STARTUP_REJECT" : "RUNTIME_STOP";
	FlushProcessingForGate();
	const auto started = std::chrono::steady_clock::now();
	const double deadline = std::chrono::duration<double>(
		std::chrono::system_clock::now().time_since_epoch()).count() + timeout;
	bool snapshot = false;
	if (!startup && m_dftSnapshots)
	{
		snapshot = WriteDFTSnapshot(id, energy, maxEnergy);
	}
	TiXmlDocument document;
	TiXmlElement* header = new TiXmlElement(startup ? "StartupGate" : "RuntimeGate");
	header->SetAttribute("id", id);
	header->SetAttribute("runID", m_runtimeGateID);
	header->SetAttribute("reason", reason);
	header->SetAttribute("thresholdDb", threshold);
	header->SetAttribute("timestep", std::to_string(FDTD_Eng->GetNumberOfTimesteps()));
	header->SetAttribute("timestepSeconds", Precise(FDTD_Op->GetTimestep()));
	header->SetAttribute("timeSeconds", Precise(FDTD_Eng->GetNumberOfTimesteps() * FDTD_Op->GetTimestep()));
	header->SetAttribute("iterationSeconds", Precise(iterationSeconds));
	header->SetAttribute("energyDb", maxEnergy > 0 && energy > 0 ? Precise(10 * std::log10(energy / maxEnergy)) : "0");
	header->SetAttribute("timeoutSeconds", std::to_string(timeout));
	header->SetAttribute("deadlineUnixSeconds", Precise(deadline));
	header->SetAttribute("maxTimesteps", std::to_string(NrTS));
	header->SetAttribute("endCriteria", Precise(endCrit));
	header->SetAttribute("snapshotAvailable", snapshot ? "true" : "false");
	header->SetAttribute("converged", "false");
	document.LinkEndChild(header);
	auto publish = [&](const char* state)
	{
		header->SetAttribute("state", state);
		if (!document.SaveFile(temporary) || std::rename(temporary, witness) != 0)
		{
			throw std::runtime_error("Run gate: cannot publish witness");
		}
	};
	if (reason == "decay" && !snapshot)
	{
		publish("snapshot-failed");
		return true;
	}
	publish("paused");
	std::cout << (startup ? "Startup gate: waiting id=" : "Runtime gate: paused id=") << id
		<< " at timestep " << FDTD_Eng->GetNumberOfTimesteps() << std::endl;
	while (true)
	{
		if (CheckAbortCond())
		{
			if (!startup)
			{
				publish("stopped");
			}
			return false;
		}
		if (m_dftSnapshots)
		{
			CheckDFTSnapshot(energy, maxEnergy);
		}
		if (Decision(stop, id))
		{
			if (startup)
			{
				throw std::runtime_error("Startup gate: qualification rejected");
			}
			publish("stopped");
			return false;
		}
		if (Decision(proceed, id))
		{
			if (!startup)
			{
				publish("resumed");
			}
			std::cout << (startup ? "Startup gate: accepted; continuing without reinitialization"
				: "Runtime gate: resumed without reinitialization") << std::endl;
			return true;
		}
		if (std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count() >= timeout)
		{
			if (startup)
			{
				throw std::runtime_error("Startup gate: qualification timed out");
			}
			publish("timed-out");
			std::cerr << "Runtime gate: controller timed out; resuming FDTD " << id << std::endl;
			return true;
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
}

bool openEMS::CheckRuntimeGates(double energy, double maxEnergy, double iterationSeconds)
{
	if (m_runtimeControl && Exists("RUNTIME_PAUSE"))
	{
		std::ifstream file("RUNTIME_PAUSE", std::ios::binary);
		std::string id;
		std::getline(file, id);
		file.close();
		const std::string prefix = m_runtimeGateID + "-";
		if (id.size() > 64 || id.compare(0, prefix.size(), prefix) != 0 ||
			id.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_") != std::string::npos)
		{
			throw std::runtime_error("Runtime gate: invalid manual pause ID");
		}
		std::remove("RUNTIME_PAUSE");
		if (!WaitRunGate(id, "manual", 0, energy, maxEnergy, iterationSeconds))
		{
			return false;
		}
	}
	if (!(maxEnergy > 0 && energy > 0) ||
		FDTD_Eng->GetNumberOfTimesteps() < m_runtimeExcitationEnd)
	{
		return true;
	}
	const double decay = 10.0 * std::log10(energy / maxEnergy);
	while (m_runtimeGateIndex < m_runtimeGateLevels.size() && decay <= m_runtimeGateLevels[m_runtimeGateIndex])
	{
		int threshold = m_runtimeGateLevels[m_runtimeGateIndex++];
		if (!WaitRunGate(m_runtimeGateID + "-d" + std::to_string(-threshold),
			"decay", threshold, energy, maxEnergy, iterationSeconds))
		{
			return false;
		}
	}
	return true;
}

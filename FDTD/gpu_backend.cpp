/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "gpu_backend.h"
#include "gpu_backend_reference.h"
#include "gpu_plugin.h"
#ifdef OPENEMS_WITH_DOMAIN_RUNTIME
#include "gpu_backend_domains.h"
#endif
#include <memory>
#include <set>
#include <limits>
#ifdef OPENEMS_WITH_HIP
#include "hip/gpu_backend_hip.h"
#endif
#ifdef OPENEMS_WITH_METAL
#include "metal/gpu_backend_metal.h"
#endif

namespace
{
struct GPU_Device
{
	const GPU_Plugin* plugin;
	int index;              //!< device number within the plugin
	GPU_DeviceInfo info;
};

//! Discrete GPUs first, then the one with more memory
bool BetterDevice(const GPU_Device& a, const GPU_Device& b)
{
	if (a.info.integrated!=b.info.integrated)
		return !a.info.integrated;
	return a.info.memory>b.info.memory;
}

std::string DeviceLine(size_t n, const GPU_Device& dev)
{
	std::ostringstream line;
	line << n << ": " << dev.info.name << " (" << dev.plugin->name << ", "
	     << std::fixed << std::setprecision(1) << dev.info.memory/1073741824.0 << " GiB"
	     << (dev.info.integrated ? ", integrated" : "") << ")";
	return line.str();
}
}

GPU_Backend* GPU_Backend::New(const std::string& name)
{
	std::vector<size_t> selected;
	std::vector<uint32_t> cuts;
	if (g_settings.hasOption("gpu-domain-cuts") && !g_settings.hasOption("gpu-devices"))
	{
		throw std::invalid_argument("--gpu-domain-cuts requires --gpu-devices");
	}
	if (g_settings.hasOption("gpu-devices"))
	{
#ifndef OPENEMS_WITH_DOMAIN_RUNTIME
		throw std::runtime_error("--gpu-devices requires ENABLE_DOMAIN_RUNTIME");
#else
		if (g_settings.GetGPUDevice() >= 0)
		{
			throw std::invalid_argument("choose either --gpu-device or --gpu-devices");
		}
		const std::string value = g_settings.getOption("gpu-devices").as<std::string>();
		std::stringstream input(value);
		std::string token;
		std::set<size_t> unique;
		while (std::getline(input, token, ','))
		{
			if (token.empty() || token.find_first_not_of("0123456789") != std::string::npos)
			{
				throw std::invalid_argument(
				    "--gpu-devices requires unique nonnegative comma-separated indices");
			}
			const size_t index = std::stoull(token);
			if (!unique.insert(index).second)
			{
				throw std::invalid_argument("duplicate --gpu-devices index");
			}
			selected.push_back(index);
		}
		if (selected.size() < 2 || value.back() == ',')
		{
			throw std::invalid_argument("--gpu-devices requires at least two indices");
		}
		if (g_settings.hasOption("gpu-domain-cuts"))
		{
			const std::string boundaries =
			    g_settings.getOption("gpu-domain-cuts").as<std::string>();
			std::stringstream list(boundaries);
			while (std::getline(list, token, ','))
			{
				if (token.empty() || token.find_first_not_of("0123456789") != std::string::npos)
				{
					throw std::invalid_argument(
					    "--gpu-domain-cuts requires nonnegative integer boundaries");
				}
				const auto cut = std::stoull(token);
				if (cut > std::numeric_limits<uint32_t>::max())
				{
					throw std::out_of_range("--gpu-domain-cuts boundary exceeds uint32 range");
				}
				cuts.push_back(static_cast<uint32_t>(cut));
			}
			if (cuts.size() != selected.size() + 1 || boundaries.back() == ',')
			{
				throw std::invalid_argument(
				    "--gpu-domain-cuts requires one more boundary than devices");
			}
		}
#endif
	}
	if (name=="reference")
	{
#ifdef OPENEMS_WITH_DOMAIN_RUNTIME
		if (!selected.empty())
		{
			std::vector<std::unique_ptr<GPU_Backend>> domains;
			for (size_t index : selected)
			{
				UNUSED(index);
				domains.emplace_back(new GPU_Backend_Reference());
			}
			return new GPU_Backend_Domains(std::move(domains), cuts);
		}
#endif
		return new GPU_Backend_Reference();
	}
	if (name!="auto")
		std::cerr << "GPU_Backend::New: Warning: unknown GPU backend \"" << name << "\", using the best available" << std::endl;

	// every backend there is, built in or a plugin; the device numbers of --gpu-device
	// follow this order
	std::vector<std::string> messages;
	std::vector<const GPU_Plugin*> plugins;
#ifdef OPENEMS_WITH_HIP
	plugins.push_back(GPU_Backend_HIP::Plugin());   // OPENEMS_GPU_PLUGINS=OFF
#endif
	const std::vector<const GPU_Plugin*> loaded = GPU_LoadPlugins(messages);
	plugins.insert(plugins.end(), loaded.begin(), loaded.end());
#ifdef OPENEMS_WITH_METAL
	plugins.push_back(GPU_Backend_Metal::Plugin());
#endif

	std::vector<GPU_Device> devices;
	for (const GPU_Plugin* plugin : plugins)
	{
		const int count = plugin->DeviceCount();
		if (count==0)
			messages.push_back(std::string("no ") + plugin->name + " device");
		for (int n=0; n<count; ++n)
		{
			GPU_Device dev;
			dev.plugin = plugin;
			dev.index = n;
			if (plugin->DeviceInfo(n, &dev.info))
				devices.push_back(dev);
		}
	}

	const int chosen = g_settings.GetGPUDevice();
	if (g_settings.GetVerboseLevel()>0)
	{
		std::cout << "GPU devices:" << std::endl;
		for (size_t n=0; n<devices.size(); ++n)
			std::cout << "  " << DeviceLine(n, devices.at(n)) << std::endl;
	}
	if ((chosen<-1) || (chosen>=(int)devices.size()))
	{
		std::ostringstream err;
		err << "--gpu-device=" << chosen << ", but there " << (devices.size()==1 ? "is " : "are ")
		    << devices.size() << " GPU device" << (devices.size()==1 ? "" : "s");
		for (size_t n=0; n<devices.size(); ++n)
			err << (n==0 ? ": " : ", ") << DeviceLine(n, devices.at(n));
		throw std::runtime_error(err.str());
	}

#ifdef OPENEMS_WITH_DOMAIN_RUNTIME
	if (!selected.empty())
	{
		std::vector<std::unique_ptr<GPU_Backend>> domains;
		for (size_t index : selected)
		{
			if (index >= devices.size())
			{
				throw std::out_of_range("--gpu-devices index exceeds available device count");
			}
			const auto& device = devices[index];
			std::unique_ptr<GPU_Backend> backend(device.plugin->NewBackend(device.index));
			if (!backend)
			{
				throw std::runtime_error("cannot initialize requested domain device");
			}
			domains.push_back(std::move(backend));
		}
		return new GPU_Backend_Domains(std::move(domains), cuts);
	}
#endif

	// the chosen device, else the best one that can run
	std::vector<size_t> order;
	if (chosen>=0)
		order.push_back(chosen);
	else
	{
		for (size_t n=0; n<devices.size(); ++n)
			order.push_back(n);
		std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {return BetterDevice(devices.at(a), devices.at(b));});
	}
	for (size_t n : order)
	{
		const GPU_Device& dev = devices.at(n);
		if (GPU_Backend* backend = dev.plugin->NewBackend(dev.index))
			return backend;
		messages.push_back("cannot run on device " + DeviceLine(n, dev));
	}

	for (const std::string& msg : messages)
		std::cerr << "GPU_Backend::New: Warning: " << msg << std::endl;
	// not the reference backend: it is far slower than the CPU engines, the caller has to say so
	return NULL;
}


#include "gpu_cell_model.h"

bool GPU_Backend::InitDomain(const Operator* op, const GPU_CellModel& model, unsigned int begin, unsigned int end)
{
	GPU_CellModel local(model, begin, end);
	m_Domain = true;
	m_Begin = begin;
	m_End = end;
	m_Plane = uint64_t(model.numLines[1]) * model.numLines[2];
	m_GlobalNodes = uint64_t(model.numLines[0]) * m_Plane;
	m_LocalNodes = uint64_t(local.numLines[0]) * m_Plane;
	return Init(op, local);
}

uint64_t GPU_Backend::LocalFlat(uint64_t global) const
{
	if (!m_Domain)
	{
		return global;
	}
	const uint64_t node = global % m_GlobalNodes;
	const unsigned int x = unsigned(node / m_Plane);
	if (!OwnsX(x))
	{
		return 0; // filtered out by GPU_PointList before device upload
	}
	return (global / m_GlobalNodes) * m_LocalNodes + (x - m_Begin + 1) * m_Plane + node % m_Plane;
}

bool GPU_Backend::DownloadRange(bool currents, size_t offset, size_t count, float* dst)
{
	DeviceDownload(dst, DeviceField(currents) + offset, count * sizeof(float));
	return true;
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "cartesian_subgrid_runtime.h"
#include "tools/global.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
std::vector<std::string> Split(const std::string &text)
{
	std::vector<std::string> result;
	if (text.empty())
	{
		return result;
	}
	std::size_t start = 0;
	while (true)
	{
		const auto end = text.find(',', start);
		const auto item = text.substr(start, end == std::string::npos ? end : end - start);
		if (item.empty() || item.find_first_of(" \t\r\n") != std::string::npos)
		{
			throw std::invalid_argument("CartesianSubgrid: empty or whitespace domain item");
		}
		result.push_back(item);
		if (end == std::string::npos)
		{
			return result;
		}
		start = end + 1;
	}
}
struct RestoreDevice
{
	int previous = g_settings.GetGPUDevice();
	~RestoreDevice() { g_settings.SetGPUDevice(previous); }
};
}

CartesianSubgridRuntime CartesianSubgridRuntime::Parse(const std::string &axis,
	const std::string &cuts, const std::string &workers)
{
	if (axis.size() != 1 || axis[0] < 'x' || axis[0] > 'z')
	{
		throw std::invalid_argument("CartesianSubgrid: domain axis must be x, y or z");
	}
	CartesianSubgridRuntime result;
	result.axis = static_cast<unsigned int>(axis[0] - 'x');
	for (const auto &text : Split(cuts == "none" ? "" : cuts))
	{
		std::size_t consumed = 0;
		const double value = std::stod(text, &consumed);
		if (consumed != text.size() || !std::isfinite(value) ||
			(!result.drawingCuts.empty() && value <= result.drawingCuts.back()))
		{
			throw std::invalid_argument("CartesianSubgrid: domain cuts must be finite and strictly increasing");
		}
		result.drawingCuts.push_back(value);
	}
	for (const auto &text : Split(workers))
	{
		const auto at = text.find('@');
		Worker worker{text.substr(0, at), -1};
		if (worker.backend != "auto" && worker.backend != "reference")
		{
			throw std::invalid_argument("CartesianSubgrid: worker must be auto[@device] or reference");
		}
		if (at != std::string::npos)
		{
			const auto number = text.substr(at + 1);
			if (worker.backend != "auto" || number.empty() ||
				!std::all_of(number.begin(), number.end(), [](char c) { return c >= '0' && c <= '9'; }))
			{
				throw std::invalid_argument("CartesianSubgrid: only auto workers accept a nonnegative device index");
			}
			const auto index = std::stoull(number);
			if (index > static_cast<unsigned int>(std::numeric_limits<int>::max()))
			{
				throw std::invalid_argument("CartesianSubgrid: worker device index overflows");
			}
			worker.device = static_cast<int>(index);
		}
		result.workers.push_back(std::move(worker));
	}
	if (result.workers.empty() || result.workers.size() > 64 ||
		result.workers.size() != result.drawingCuts.size() + 1)
	{
		throw std::invalid_argument("CartesianSubgrid: one through 64 workers require exactly one fewer domain cuts");
	}
	return result;
}

boost::program_options::options_description CartesianSubgridRuntime::Options()
{
	namespace po = boost::program_options;
	po::options_description options("Experimental Cartesian-subgrid domains");
	options.add_options()
		("subgrid-domain-axis", po::value<std::string>(), "Fixed assignment axis: x, y or z")
		("subgrid-domain-cuts", po::value<std::string>(), "Comma-separated global drawing-coordinate cuts; none for one worker")
		("subgrid-domain-workers", po::value<std::string>(), "Comma-separated auto[@device] or explicit reference workers");
	return options;
}

std::optional<CartesianSubgridRuntime> CartesianSubgridRuntime::ReadOptions()
{
	const bool axis = g_settings.hasOption("subgrid-domain-axis");
	const bool cuts = g_settings.hasOption("subgrid-domain-cuts");
	const bool workers = g_settings.hasOption("subgrid-domain-workers");
	if (!axis && !cuts && !workers)
	{
		return std::nullopt;
	}
	if (!axis || !cuts || !workers)
	{
		throw std::invalid_argument("CartesianSubgrid: provide axis, cuts and workers together");
	}
	return Parse(g_settings.getOption("subgrid-domain-axis").as<std::string>(),
		g_settings.getOption("subgrid-domain-cuts").as<std::string>(),
		g_settings.getOption("subgrid-domain-workers").as<std::string>());
}

std::vector<std::unique_ptr<GPU_Backend>> CartesianSubgridRuntime::CreateBackends(bool referenceOnly) const
{
	RestoreDevice restore;
	std::vector<std::unique_ptr<GPU_Backend>> result;
	for (const auto &worker : workers)
	{
		if (referenceOnly && worker.backend != "reference")
		{
			throw std::invalid_argument("CartesianSubgrid: gpu-reference requires explicitly reference workers");
		}
		g_settings.SetGPUDevice(worker.device < 0 ? restore.previous : worker.device);
		result.emplace_back(GPU_Backend::New(worker.backend));
		if (!result.back())
		{
			throw std::invalid_argument("CartesianSubgrid: requested worker backend is unavailable");
		}
	}
	return result;
}

/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef CARTESIAN_SUBGRID_RUNTIME_H
#define CARTESIAN_SUBGRID_RUNTIME_H

#include "gpu_backend.h"
#include <boost/program_options.hpp>
#include <memory>
#include <optional>
#include <string>
#include <vector>

//! Execution placement only. Geometry XML and the external partition planner
//! remain independent of physical devices and transport implementation.
struct CartesianSubgridRuntime
{
	struct Worker
	{
		std::string backend;
		int device = -1;
	};
	unsigned int axis = 0;
	std::vector<double> drawingCuts;
	std::vector<Worker> workers;
	static CartesianSubgridRuntime Parse(const std::string &axis,
		const std::string &cuts, const std::string &workers);
	static boost::program_options::options_description Options();
	static std::optional<CartesianSubgridRuntime> ReadOptions();
	std::vector<std::unique_ptr<GPU_Backend>> CreateBackends(bool referenceOnly) const;
};
#endif

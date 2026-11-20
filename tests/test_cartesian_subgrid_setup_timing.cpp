/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "FDTD/cartesian_subgrid.h"
#include "FDTD/subgrid_setup_timing.h"
#include "tools/constants.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

static void Require(bool value, const char *message)
{
	if (!value)
	{
		throw std::runtime_error(message);
	}
}

class ProfileEnvironment
{
  public:
	explicit ProfileEnvironment(const char *value)
	{
		const char *old = std::getenv("OPENEMS_SUBGRID_SETUP_PROFILE");
		m_Present = old != nullptr;
		m_Previous = old ? old : "";
		Set(value);
	}
	~ProfileEnvironment()
	{
		Set(m_Present ? m_Previous.c_str() : nullptr);
	}
  private:
	static void Set(const char *value)
	{
#ifdef _WIN32
		_putenv_s("OPENEMS_SUBGRID_SETUP_PROFILE", value ? value : "");
#else
		if (value)
		{
			setenv("OPENEMS_SUBGRID_SETUP_PROFILE", value, 1);
		}
		else
		{
			unsetenv("OPENEMS_SUBGRID_SETUP_PROFILE");
		}
#endif
	}
	bool m_Present;
	std::string m_Previous;
};

class CaptureLog
{
  public:
	CaptureLog() : m_Previous(std::clog.rdbuf(m_Output.rdbuf())) {}
	~CaptureLog()
	{
		std::clog.rdbuf(m_Previous);
	}
	std::string Text() const
	{
		return m_Output.str();
	}
  private:
	std::ostringstream m_Output;
	std::streambuf *m_Previous;
};

struct State
{
	CartesianSubgridDeviceData data;
	double step = 0, energy = 0;
	std::string log;
};

static State Run(const char *value)
{
	ProfileEnvironment environment(value);
	CaptureLog log;
	CartesianSubgrid grid({6, 6, 6}, {.001, .001, .001}, {{2, 2, 2}, {4, 4, 4}}, EPS0, MUE0);
	grid.SetElectric([](unsigned int n, const CartesianSubgrid::Point &p)
	{
		return (n + 1) * std::sin(100 * p[0]) * std::cos(70 * p[1]) * std::sin(90 * p[2]);
	});
	grid.SetMagnetic([](unsigned int n, const CartesianSubgrid::Point &p)
	{
		return (n + 1) * .001 * std::cos(50 * p[0]) * std::sin(80 * p[1]) * std::cos(60 * p[2]);
	});
	const double step = .5 * grid.StableTimestep();
	for (unsigned int i = 0; i < 64; ++i)
	{
		grid.Step(step);
	}
	return {grid.DeviceStateData(), grid.StableTimestep(), grid.Energy(), log.Text()};
}

int main()
{
	try
	{
		const auto absent = Run(nullptr), disabled = Run("0"), invalid = Run("yes"), enabled = Run("1");
		Require(absent.log.empty() && disabled.log.empty() && invalid.log.empty(), "setup timing is not opt-in");
		for (const State *state : {&disabled, &invalid, &enabled})
		{
			Require(absent.data.electric == state->data.electric && absent.data.magnetic == state->data.magnetic,
				"setup profiling changed propagated fields");
			Require(absent.step == state->step && absent.energy == state->energy,
				"setup profiling changed stability or energy");
		}
		const char *names[] = {"leaf_addresses", "leaf_materials_and_opaque", "opaque_edge_elimination",
			"pec_classification", "master_edges_and_constraints", "electric_mass_assembly", "interface_factor", "curl_and_step_bound"};
		for (const char *name : names)
		{
			Require(enabled.log.find(std::string("group=grid name=") + name + " seconds=") != std::string::npos,
				"missing completed preparation phase");
		}
		Require(enabled.log.find("completed_stages=8 complete=1") != std::string::npos, "missing successful preparation marker");
		{
			ProfileEnvironment environment("1");
			CaptureLog log;
			CartesianSubgrid::Axes axes;
			for (auto &axis : axes)
			{
				for (unsigned int i = 0; i <= 6; ++i)
				{
					axis.push_back(.001 * i);
				}
			}
			bool rejected = false;
			try
			{
				CartesianSubgrid invalidGrid(axes, EPS0, MUE0,
					{{{{2, 2, 2}, {4, 4, 4}}, {2, 2, 2}}}, false,
					[](const CartesianSubgrid::Point &, const CartesianSubgrid::Point &) -> CartesianSubgrid::Medium
					{
						throw std::invalid_argument("controlled material failure");
					});
			}
			catch (const std::invalid_argument &error)
			{
				rejected = std::string(error.what()) == "controlled material failure";
			}
			Require(rejected, "timing replaced the original preparation failure");
			Require(log.Text().find("completed_stages=1 complete=0") != std::string::npos,
				"failed preparation was reported as complete");
		}
		{
			std::ostringstream failed;
			failed.exceptions(std::ios::badbit | std::ios::failbit);
			try
			{
				failed.setstate(std::ios::badbit);
			}
			catch (const std::ios_base::failure &)
			{
			}
			CartesianSubgridSetupTiming timing("failed_stream", true, failed);
			timing.Phase("diagnostic_failure");
			timing.Complete();
		}
		std::cout << "SETUP_TIMING phases=8 exact_fields=1 exact_step_energy=1 incomplete_guard=1 stream_failure=1 opt_in=1\n";
		std::cout << "PASS: setup timing preserves fields, disabled output and numerical exceptions\n";
		return 0;
	}
	catch (const std::exception &error)
	{
		std::cerr << "FAIL: " << error.what() << '\n';
		return 1;
	}
}

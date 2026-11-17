/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef SUBGRID_SETUP_TIMING_H
#define SUBGRID_SETUP_TIMING_H

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>

//! Optional serial preparation diagnostics, never called by a field update.
class CartesianSubgridSetupTiming
{
  public:
	static bool Requested()
	{
		const char *value = std::getenv("OPENEMS_SUBGRID_SETUP_PROFILE");
		return value && std::strcmp(value, "1") == 0;
	}

	explicit CartesianSubgridSetupTiming(const char *group, bool enabled = Requested(),
		std::ostream &output = std::clog)
		: m_Group(group), m_Enabled(enabled), m_Output(output)
	{
		if (m_Enabled)
		{
			m_Start = m_Last = Clock::now();
		}
	}

	~CartesianSubgridSetupTiming() noexcept
	{
		if (m_Enabled)
		{
			try
			{
				m_Output << "CartesianSubgrid setup total group=" << m_Group
					<< " seconds=" << Seconds(Clock::now() - m_Start)
					<< " completed_stages=" << m_Stages << " complete=" << m_Complete << '\n';
			}
			catch (...)
			{
				// A diagnostic stream must not replace a numerical exception.
			}
		}
	}

	void Phase(const char *name) noexcept
	{
		if (m_Enabled)
		{
			const auto now = Clock::now();
			try
			{
				m_Output << "CartesianSubgrid setup phase group=" << m_Group << " name=" << name
					<< " seconds=" << Seconds(now - m_Last) << '\n';
			}
			catch (...)
			{
				// Keep diagnostic failures separate from solver behavior.
			}
			m_Last = now;
			++m_Stages;
		}
	}

	void Complete() noexcept
	{
		m_Complete = true;
	}

  private:
	using Clock = std::chrono::steady_clock;
	static double Seconds(Clock::duration duration)
	{
		return std::chrono::duration<double>(duration).count();
	}
	const char *m_Group;
	bool m_Enabled, m_Complete = false;
	std::ostream &m_Output;
	Clock::time_point m_Start{}, m_Last{};
	unsigned int m_Stages = 0;
};
#endif

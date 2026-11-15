// Optional, coarse host-wall timings of setup, never of the FDTD update loop.
// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_SETUP_PROFILE_H
#define OPENEMS_SETUP_PROFILE_H

#include "openems_global.h"
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

class OPENEMS_EXPORT SetupProfile
{
public:
	using Clock = double (*)();
	struct Entry
	{
		std::string name;
		std::string detail;
		int parent;
		double seconds;
	};

	// The context is local to the calling setup thread. Nested setup sessions
	// restore their predecessor; worker threads do not time individual cells.
	SetupProfile(bool enabled, std::ostream& output, Clock clock = SteadySeconds);
	~SetupProfile() noexcept;
	SetupProfile(const SetupProfile&) = delete;
	SetupProfile& operator=(const SetupProfile&) = delete;
	void Complete(bool operatorOnly = false);
	static bool Active();
	static bool FromEnvironment();
	static double SteadySeconds();

	class OPENEMS_EXPORT Phase
	{
	public:
		explicit Phase(const char* name, const std::string& detail = std::string());
		~Phase() noexcept;
		Phase(const Phase&) = delete;
		Phase& operator=(const Phase&) = delete;
	private:
		SetupProfile* m_Owner;
		size_t m_Index;
		double m_Start;
	};

private:
	static thread_local SetupProfile* s_Current;
	bool m_Enabled;
	bool m_Complete;
	bool m_OperatorOnly;
	SetupProfile* m_Previous;
	std::ostream& m_Output;
	Clock m_Clock;
	double m_Start;
	std::vector<Entry> m_Entries;
	std::vector<size_t> m_Stack;
	void Write(double seconds);
};
#endif

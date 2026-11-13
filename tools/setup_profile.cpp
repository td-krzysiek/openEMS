// SPDX-License-Identifier: GPL-3.0-or-later
#include "setup_profile.h"
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <locale>
#include <ostream>
#include <sstream>

thread_local SetupProfile* SetupProfile::s_Current = nullptr;

namespace
{
void Quoted(std::ostream& out, const std::string& text)
{
	out << '"';
	for (unsigned char c : text)
	{
		if (c == '"' || c == '\\')
		{
			out << '\\' << char(c);
		}
		else if (c < 0x20)
		{
			const char* hex = "0123456789abcdef";
			out << "\\u00" << hex[c >> 4] << hex[c & 15];
		}
		else
		{
			out << char(c);
		}
	}
	out << '"';
}
}

SetupProfile::SetupProfile(bool enabled, std::ostream& output, Clock clock)
	: m_Enabled(enabled), m_Complete(false), m_OperatorOnly(false),
	  m_Previous(s_Current), m_Output(output), m_Clock(clock),
	  m_Start(enabled ? clock() : 0)
{
	s_Current = enabled ? this : nullptr;
}

SetupProfile::~SetupProfile() noexcept
{
	s_Current = m_Previous;
	if (m_Enabled)
	{
		try
		{
			Write(m_Clock() - m_Start);
		}
		catch (...)
		{
			// A diagnostic output failure must not terminate a solver on exit.
		}
	}
}

bool SetupProfile::Active()
{
	return s_Current != nullptr;
}

bool SetupProfile::FromEnvironment()
{
	const char* value = std::getenv("OPENEMS_PROFILE_SETUP");
	return value && std::strcmp(value, "1") == 0;
}

double SetupProfile::SteadySeconds()
{
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

void SetupProfile::Complete(bool operatorOnly)
{
	m_Complete = true;
	m_OperatorOnly = operatorOnly;
}

SetupProfile::Phase::Phase(const char* name, const std::string& detail)
	: m_Owner(s_Current), m_Index(0), m_Start(0)
{
	if (m_Owner)
	{
		m_Index = m_Owner->m_Entries.size();
		m_Owner->m_Entries.push_back({name, detail,
			m_Owner->m_Stack.empty() ? -1 : int(m_Owner->m_Stack.back()), 0});
		m_Owner->m_Stack.push_back(m_Index);
		m_Start = m_Owner->m_Clock();
	}
}

SetupProfile::Phase::~Phase() noexcept
{
	if (m_Owner)
	{
		m_Owner->m_Entries[m_Index].seconds = m_Owner->m_Clock() - m_Start;
		m_Owner->m_Stack.pop_back();
	}
}

void SetupProfile::Write(double seconds)
{
	// Build one record without changing the formatting/locale of the solver log.
	std::ostringstream out;
	out.imbue(std::locale::classic());
	out << std::setprecision(17)
	    << "Setup profile: {\"schema\":\"openems.setup-profile.v1\",\"status\":\""
	    << (m_Complete ? "complete" : "incomplete") << "\",\"scope\":\""
	    << (m_OperatorOnly ? "operator-only" : "setup")
	    << "\",\"seconds\":" << seconds << ",\"phases\":[";
	for (size_t i = 0; i < m_Entries.size(); ++i)
	{
		const Entry& e = m_Entries[i];
		if (i)
		{
			out << ',';
		}
		out << "{\"name\":";
		Quoted(out, e.name);
		out << ",\"detail\":";
		Quoted(out, e.detail);
		out << ",\"parent\":" << e.parent << ",\"seconds\":" << e.seconds << '}';
	}
	out << "]}\n";
	m_Output << out.str();
}

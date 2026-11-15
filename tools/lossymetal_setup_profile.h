// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LOSSYMETAL_SETUP_PROFILE_H
#define LOSSYMETAL_SETUP_PROFILE_H

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

// Setup only. Samples are inclusive, not an additive phase decomposition.
class LossyMetalSetupProfile
{
public:
	using Clock = double (*)();
	enum Query {ColumnCandidates, EdgeProperty, Contributions, IsMetal, IsLossy,
		InsideLength, Extent, LineCandidates, FitLookup, FitSolve, QueryCount};
	struct Stats
	{
		uint64_t calls = 0;
		uint64_t samples = 0;
		double seconds = 0;
	};
	struct Phase
	{
		std::string name;
		double seconds;
	};

	static double Now()
	{
		return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
	}

	explicit LossyMetalSetupProfile(bool enabled, std::ostream& output = std::cerr,
		uint64_t interval = 1024, Clock clock = Now)
		: m_Enabled(enabled), m_Output(output), m_Interval(interval ? interval : 1024),
		  m_Clock(clock), m_Previous(s_Current), m_Start(enabled ? clock() : 0), m_Last(m_Start)
	{
		s_Current = enabled ? this : nullptr;
	}

	~LossyMetalSetupProfile() noexcept
	{
		s_Current = m_Previous;
		if (m_Enabled) {
			try {
				Write();
			} catch (...) {
				// Optional diagnostics must not fail the numerical setup.
			}
		}
	}

	LossyMetalSetupProfile(const LossyMetalSetupProfile&) = delete;
	LossyMetalSetupProfile& operator=(const LossyMetalSetupProfile&) = delete;

	static bool FromEnvironment()
	{
		const char* value = std::getenv("OPENEMS_PROFILE_LOSSY_SETUP");
		return value && std::strcmp(value, "1") == 0;
	}

	void FinishPhase(const char* name)
	{
		if (m_Enabled) {
			const double now = m_Clock();
			m_Phases.push_back({name, now - m_Last});
			m_Last = now;
		}
	}

	void Complete() {m_Complete = true;}
	const Stats& GetStats(Query query) const {return m_Stats[query];}

	class Sample
	{
	public:
		explicit Sample(Query query) : m_Owner(s_Current), m_Query(query)
		{
			if (m_Owner) {
				Stats& stats = m_Owner->m_Stats[query];
				++stats.calls;
				// Fit misses are rare and expensive: measure every solve.
				m_Sampled = query == FitSolve || (stats.calls - 1) % m_Owner->m_Interval == 0;
				if (m_Sampled) {
					m_Start = m_Owner->m_Clock();
				}
			}
		}

		~Sample()
		{
			if (m_Sampled) {
				Stats& stats = m_Owner->m_Stats[m_Query];
				stats.seconds += m_Owner->m_Clock() - m_Start;
				++stats.samples;
			}
		}

		Sample(const Sample&) = delete;
		Sample& operator=(const Sample&) = delete;

	private:
		LossyMetalSetupProfile* m_Owner;
		Query m_Query;
		bool m_Sampled = false;
		double m_Start = 0;
	};

private:
	inline static thread_local LossyMetalSetupProfile* s_Current = nullptr;
	bool m_Enabled;
	bool m_Complete = false;
	std::ostream& m_Output;
	uint64_t m_Interval;
	Clock m_Clock;
	LossyMetalSetupProfile* m_Previous;
	double m_Start, m_Last;
	std::array<Stats, QueryCount> m_Stats;
	std::vector<Phase> m_Phases;

	void Write() const
	{
		static const char* names[] = {"column_candidates", "edge_property", "surface_contributions",
			"is_metal", "is_lossy", "inside_length", "extent", "line_candidates", "fit_lookup", "fit_solve"};
		std::ostringstream out;
		out.imbue(std::locale::classic());
		out << std::setprecision(17) << "LossyMetal setup profile: {\"schema\":\"openems.lossymetal-setup.v1\","
			<< "\"status\":\"" << (m_Complete ? "complete" : "incomplete") << "\",\"seconds\":"
			<< m_Clock() - m_Start << ",\"sample_interval\":" << m_Interval << ",\"phases\":[";
		for (size_t i = 0; i < m_Phases.size(); ++i) {
			if (i) {
				out << ',';
			}
			out << "{\"name\":\"" << m_Phases[i].name << "\",\"seconds\":" << m_Phases[i].seconds << '}';
		}
		out << "],\"queries\":[";
		for (size_t i = 0; i < QueryCount; ++i) {
			const Stats& stats = m_Stats[i];
			if (i) {
				out << ',';
			}
			out << "{\"name\":\"" << names[i] << "\",\"calls\":" << stats.calls
				<< ",\"samples\":" << stats.samples << ",\"sampled_seconds\":" << stats.seconds
				<< ",\"estimated_inclusive_seconds\":" << (stats.samples ? stats.seconds * double(stats.calls) / stats.samples : 0) << '}';
		}
		out << "]}\n";
		m_Output << out.str();
	}
};
#endif

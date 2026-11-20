// SPDX-License-Identifier: GPL-3.0-or-later
#include "tools/lossymetal_setup_profile.h"
#include <sstream>
#include <thread>

namespace
{
int failures = 0;
int clocks = 0;
double FakeClock()
{
	return clocks++;
}
void Check(bool ok, const char* message)
{
	if (!ok) {
		std::cerr << "FAIL: " << message << '\n';
		++failures;
	}
}
}

int main()
{
	using Profile = LossyMetalSetupProfile;
	std::ostringstream disabled;
	{
		Profile profile(false, disabled, 4, FakeClock);
		Profile::Sample query(Profile::IsLossy);
		profile.FinishPhase("off");
	}
	Check(clocks == 0 && disabled.str().empty(), "disabled has no clock reads or output");
	std::ostringstream out;
	{
		Profile profile(true, out, 4, FakeClock);
		for (int i = 0; i < 9; ++i) {
			Profile::Sample query(Profile::IsLossy);
		}
		Check(profile.GetStats(Profile::IsLossy).calls == 9, "exact call count");
		Check(profile.GetStats(Profile::IsLossy).samples == 3, "first call and every fourth call sampled");
		Check(profile.GetStats(Profile::IsLossy).seconds == 3, "sampled durations");
		{
			Profile nested(false, disabled, 4, FakeClock);
			Profile::Sample query(Profile::IsLossy);
		}
		Check(profile.GetStats(Profile::IsLossy).calls == 9, "nested disabled scope mutes queries");
		{
			Profile::Sample query(Profile::IsLossy);
		}
		Check(profile.GetStats(Profile::IsLossy).calls == 10, "outer context restored");
		std::thread worker([]() {Profile::Sample query(Profile::IsLossy);});
		worker.join();
		Check(profile.GetStats(Profile::IsLossy).calls == 10, "independent thread does not share counters");
		for (int i = 0; i < 3; ++i) {
			Profile::Sample query(Profile::FitSolve);
		}
		Check(profile.GetStats(Profile::FitSolve).samples == 3, "rare fit solves always measured");
		profile.FinishPhase("surface_scan");
		profile.Complete();
	}
	Check(out.str().find("\"status\":\"complete\"") != std::string::npos, "complete JSON record");
	Check(out.str().find("\"sample_interval\":4") != std::string::npos, "actual sampling interval recorded");
	std::ostringstream incomplete;
	{
		Profile profile(true, incomplete, 0, FakeClock);
	}
	Check(incomplete.str().find("\"status\":\"incomplete\"") != std::string::npos, "unfinished setup explicit");
	Check(incomplete.str().find("\"sample_interval\":1024") != std::string::npos, "zero interval guarded");
	{
		Profile combined(false, disabled);
		std::array<Profile::Stats, Profile::QueryCount> worker{};
		worker[Profile::Extent] = {9, 3, 2.5};
		combined.Absorb(worker);
		combined.Absorb(worker);
		const auto snapshot = combined.Snapshot();
		Check(snapshot[Profile::Extent].calls == 18 && snapshot[Profile::Extent].samples == 6 &&
			snapshot[Profile::Extent].seconds == 5, "worker counters merge without clocks or output");
	}
	return failures ? 1 : 0;
}

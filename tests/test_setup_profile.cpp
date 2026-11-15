// SPDX-License-Identifier: GPL-3.0-or-later
#include "tools/setup_profile.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace
{
double now = 0;
int calls = 0;
int failures = 0;
double FakeClock()
{
	++calls;
	return now;
}
void Check(bool condition, const char* message)
{
	if (!condition)
	{
		std::cerr << "FAIL: " << message << '\n';
		++failures;
	}
}
}

int main()
{
	std::ostringstream disabled;
	{
		SetupProfile session(false, disabled, FakeClock);
		SetupProfile::Phase phase("off");
	}
	Check(calls == 0 && disabled.str().empty(), "disabled path has no clock reads or output");
	std::ostringstream out;
	{
		SetupProfile session(true, out, FakeClock);
		{
			SetupProfile::Phase outer("outer");
			now = 2;
			{
				SetupProfile::Phase inner("inner", "quote\" slash\\ newline\n");
				now = 5;
			}
			now = 7;
		}
		std::ostringstream muted;
		{
			SetupProfile session(false, muted, FakeClock);
			Check(!SetupProfile::Active(), "disabled nested session mutes context");
		}
		Check(SetupProfile::Active(), "nested session restores context");
		bool independent = false;
		std::thread worker([&independent]() {
			independent = !SetupProfile::Active();
		});
		worker.join();
		Check(independent, "context is thread-local");
		session.Complete(true);
		now = 9;
	}
	Check(!SetupProfile::Active(), "session restores empty context");
	const std::string text = out.str();
	Check(text.find("\"status\":\"complete\",\"scope\":\"operator-only\",\"seconds\":9") != std::string::npos,
		"session status and exact elapsed time");
	Check(text.find("\"parent\":-1,\"seconds\":7") != std::string::npos, "inclusive outer duration");
	Check(text.find("\"parent\":0,\"seconds\":3") != std::string::npos, "inner parent and duration");
	Check(text.find("quote\\\" slash\\\\ newline\\u000a") != std::string::npos, "JSON string escaping");
	std::ostringstream failed;
	try
	{
		SetupProfile session(true, failed, FakeClock);
		SetupProfile::Phase phase("unwinding");
		throw std::runtime_error("test");
	}
	catch (const std::runtime_error&)
	{
	}
	Check(failed.str().find("\"status\":\"incomplete\"") != std::string::npos, "exception retains incomplete record");
	std::ostringstream broken;
	broken.exceptions(std::ios::badbit);
	try
	{
		broken.setstate(std::ios::badbit);
	}
	catch (const std::ios_base::failure&)
	{
	}
	{
		SetupProfile session(true, broken, FakeClock);
	}
	Check(!SetupProfile::Active(), "failed diagnostic output cannot retain context");
	return failures == 0 ? 0 : 1;
}

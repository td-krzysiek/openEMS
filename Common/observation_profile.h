// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_OBSERVATION_PROFILE_H
#define OPENEMS_OBSERVATION_PROFILE_H

#include <chrono>
#include <cstdint>
#include <string>
#include "openems_global.h"

//! Opt-in host-wall observation records; never inserts a device synchronization.
namespace ObservationProfile
{
OPENEMS_EXPORT bool Enabled();
OPENEMS_EXPORT void Record(const char* stage, const std::string& owner,
    unsigned int step, double seconds, std::uint64_t bytes,
    unsigned int interval = 0, unsigned int frequencies = 0);
//! Flush buffered records at processing boundaries, not in the FDTD update loop.
OPENEMS_EXPORT void Flush();

class OPENEMS_EXPORT Scope
{
public:
    Scope(const char* stage, const std::string& owner, unsigned int step,
        std::uint64_t bytes = 0, unsigned int interval = 0,
        unsigned int frequencies = 0, bool active = true);
    ~Scope();
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

private:
    bool m_Enabled;
    const char* m_Stage;
    std::string m_Owner;
    unsigned int m_Step, m_Interval, m_Frequencies;
    std::uint64_t m_Bytes;
    std::chrono::steady_clock::time_point m_Start;
};
}
#endif

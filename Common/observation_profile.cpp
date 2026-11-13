// SPDX-License-Identifier: GPL-3.0-or-later
#include "observation_profile.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>

namespace
{
void JsonString(std::ostream& stream, const std::string& value)
{
    const char hex[] = "0123456789abcdef";
    stream << '"';
    for (unsigned char c : value)
    {
        if (c == '"' || c == '\\')
        {
            stream << '\\' << c;
        }
        else if (c < 32)
        {
            stream << "\\u00" << hex[c >> 4] << hex[c & 15];
        }
        else
        {
            stream << c;
        }
    }
    stream << '"';
}

class Writer
{
public:
    Writer()
    {
        const char* path = std::getenv("OPENEMS_OBSERVATION_PROFILE");
        if (!path || !*path)
        {
            return;
        }
        m_Stream.open(path, std::ios::out | std::ios::app);
        if (!m_Stream)
        {
            std::cerr << "Observation profile: cannot open " << path << "; monitoring disabled\n";
            return;
        }
        m_Enabled = true;
        m_Stream << std::setprecision(17);
        m_Stream << "{\"type\":\"header\",\"schema\":1,\"clock\":\"host_steady_wall\","
            "\"device_waits_added\":false,\"scopes_may_overlap\":true}\n";
        m_Stream.flush();
    }

    ~Writer()
    {
        Flush();
    }

    bool Enabled() const
    {
        return m_Enabled;
    }

    void Record(const char* stage, const std::string& owner, unsigned int step,
        double seconds, std::uint64_t bytes, unsigned int interval, unsigned int frequencies)
    {
        if (!m_Enabled || !std::isfinite(seconds) || seconds < 0)
        {
            return;
        }
        std::lock_guard<std::mutex> lock(m_Mutex);
        m_Stream << "{\"type\":\"scope\",\"stage\":";
        JsonString(m_Stream, stage);
        m_Stream << ",\"owner\":";
        JsonString(m_Stream, owner);
        m_Stream << ",\"step\":" << step << ",\"host_wall_seconds\":" << seconds
            << ",\"logical_bytes\":" << bytes << ",\"sample_interval\":" << interval
            << ",\"frequencies\":" << frequencies << "}\n";
    }

    void Flush()
    {
        if (!m_Enabled)
        {
            return;
        }
        std::lock_guard<std::mutex> lock(m_Mutex);
        m_Stream.flush();
        if (!m_Stream && !m_Warned)
        {
            // A full or removed output device must not abort a simulation from a scope destructor.
            std::cerr << "Observation profile: output failed; subsequent records may be missing\n";
            m_Warned = true;
        }
    }

private:
    bool m_Enabled = false, m_Warned = false;
    std::ofstream m_Stream;
    std::mutex m_Mutex;
};

Writer& GetWriter()
{
    static Writer writer;
    return writer;
}
}

bool ObservationProfile::Enabled()
{
    return GetWriter().Enabled();
}

void ObservationProfile::Record(const char* stage, const std::string& owner,
    unsigned int step, double seconds, std::uint64_t bytes,
    unsigned int interval, unsigned int frequencies)
{
    GetWriter().Record(stage, owner, step, seconds, bytes, interval, frequencies);
}

void ObservationProfile::Flush()
{
    GetWriter().Flush();
}

ObservationProfile::Scope::Scope(const char* stage, const std::string& owner,
    unsigned int step, std::uint64_t bytes, unsigned int interval,
    unsigned int frequencies, bool active)
    : m_Enabled(active && Enabled()), m_Stage(stage), m_Step(step),
      m_Interval(interval), m_Frequencies(frequencies), m_Bytes(bytes)
{
    if (m_Enabled)
    {
        m_Owner = owner;
        m_Start = std::chrono::steady_clock::now();
    }
}

ObservationProfile::Scope::~Scope()
{
    if (m_Enabled)
    {
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - m_Start).count();
        Record(m_Stage, m_Owner, m_Step, seconds, m_Bytes, m_Interval, m_Frequencies);
    }
}

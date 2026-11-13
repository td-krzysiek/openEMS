// SPDX-License-Identifier: GPL-3.0-or-later
#include "Common/observation_profile.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <thread>
#include <vector>

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        return 2;
    }
    const std::string path = argv[1], mode = argv[2];
    std::remove(path.c_str());
    const std::string value = mode == "on" ? path : "";
#ifdef _WIN32
    _putenv_s("OPENEMS_OBSERVATION_PROFILE", value.c_str());
#else
    setenv("OPENEMS_OBSERVATION_PROFILE", value.c_str(), 1);
#endif
    if (ObservationProfile::Enabled() != (mode == "on"))
    {
        return 3;
    }
    std::vector<std::thread> workers;
    for (unsigned int t = 0; t < 4; ++t)
    {
        workers.emplace_back([t]()
        {
            for (unsigned int i = 0; i < 100; ++i)
            {
                ObservationProfile::Record("test", "quote\"slash\\newline\n", t * 100 + i,
                    .25, 1024, 7, 5);
            }
        });
    }
    for (auto& worker : workers)
    {
        worker.join();
    }
    ObservationProfile::Record("invalid", "ignored", 0, -1, 0);
    {
        ObservationProfile::Scope scope("timed", "test", 400, 4, 7, 5);
    }
    ObservationProfile::Flush();
    std::ifstream stream(path);
    if (mode != "on")
    {
        return stream.is_open() ? 4 : 0;
    }
    const std::string text((std::istreambuf_iterator<char>(stream)), {});
    size_t count = 0, pos = 0;
    while ((pos = text.find("\"type\":\"scope\"", pos)) != std::string::npos)
    {
        ++count;
        ++pos;
    }
    if (count != 401 || text.find("\\u000a") == std::string::npos ||
        text.find("quote\\\"slash\\\\newline") == std::string::npos ||
        text.find("\"stage\":\"invalid\"") != std::string::npos)
    {
        std::cerr << "Invalid observation profile count/escaping\n";
        return 5;
    }
    std::cout << "PASS optional concurrent observation records\n";
    return 0;
}

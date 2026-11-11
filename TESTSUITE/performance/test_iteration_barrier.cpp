// SPDX-License-Identifier: GPL-3.0-or-later
// Reuse and memory visibility across a fixed worker-team barrier.
#ifdef NDEBUG
#error "Compile with assertions enabled"
#endif
#include "FDTD/iteration_barrier.h"
#include <cassert>
#include <iostream>
#include <thread>
#include <vector>

int main()
{
    for (unsigned int count : {1u, 2u, 4u, 16u})
    {
        IterationBarrier barrier(count);
        std::vector<unsigned int> values(count, 0);
        std::vector<std::thread> workers;
        for (unsigned int worker=0; worker<count; ++worker)
            workers.emplace_back([&, worker]() {
                for (unsigned int round=1; round<=10000; ++round)
                {
                    values[worker] = round;
                    barrier.wait();
                    for (unsigned int value : values) assert(value==round);
                    barrier.wait(); // all reads finish before the next writes
                }
            });
        for (auto& worker : workers) worker.join();
    }
    std::cout << "Iteration barrier reuse and visibility checks passed\n";
}

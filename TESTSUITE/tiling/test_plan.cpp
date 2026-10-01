// Exact spatial/time coverage, including clipped one-half-step edge tiles.
// Compile from the repository root with:
// c++ -std=c++11 -I. TESTSUITE/tiling/test_plan.cpp tools/tiling.cpp -o /tmp/tiling-plan
#include "tools/tiling.h"
#include <iostream>
#include <stdexcept>

int main()
{
    unsigned int cases = 0;
    for (int nx=2; nx<=39; ++nx)
        for (int ny : {2, 11, 23})
            for (int size : {1, 2, 5, 10, 100})
                for (int steps : {1, 2, 5, 8})
                {
                    auto x = computeDiamondTiles1D(nx, size, 2*steps);
                    auto y = computeDiamondTiles1D(ny, size, 2*steps);
                    Tiles z;
                    z.type = TILES_RECTANGULAR;
                    z.phases = 1;
                    z.array = {{Block(2*steps, Range(0, 3))}};
                    auto plan = combineTilesTo3D(x, y, z, 2*steps, 4);
                    std::vector<int> hits(nx*ny*steps*2, 0);
                    for (const auto& team : plan)
                        for (const auto& phase : team)
                            for (const auto& tile : phase)
                                for (int eh=0; eh<2; ++eh)
                                {
                                    if (!(eh ? tile.HasCurrent() : tile.HasVoltage())) continue;
                                    const int* start = eh ? tile.currentStart : tile.voltageStart;
                                    const int* stop = eh ? tile.currentStop : tile.voltageStop;
                                    if (start[0]<0 || stop[0]>=nx || start[1]<0 || stop[1]>=ny ||
                                        start[2]!=0 || stop[2]!=3 || tile.timestep<0 || tile.timestep>=steps)
                                        throw std::runtime_error("Invalid tile bounds");
                                    for (int i=start[0]; i<=stop[0]; ++i)
                                        for (int j=start[1]; j<=stop[1]; ++j)
                                            ++hits[(2*tile.timestep+eh)*nx*ny + i*ny+j];
                                }
                    for (int hit : hits)
                        if (hit!=1)
                        {
                            std::cerr << "Invalid coverage: " << nx << ',' << ny
                                      << " tile=" << size << " steps=" << steps << '\n';
                            return 1;
                        }
                    ++cases;
                }
    std::cout << cases << " plans have exact E/H coverage\n";
}

/*
 * Copyright (C) 2023 Yifeng Li <tomli@tomli.me>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef OPENEMS_TILING_H
#define OPENEMS_TILING_H

#include <utility>
#include <vector>
#include <cstddef>

using Range = std::pair<int, int>;
using Block = std::vector<Range>;
enum { TILES_RECTANGULAR, TILES_DIAMOND };

struct Tiles
{
    int type;
    int phases;
    // Phase, spatial block, half timestep (E followed by H).
    std::vector<std::vector<Block>> array;
};

struct Range3D
{
    int timestep;
    int voltageStart[3], voltageStop[3];
    int currentStart[3], currentStop[3];
    bool HasVoltage() const { return Nonempty(voltageStart, voltageStop); }
    bool HasCurrent() const { return Nonempty(currentStart, currentStop); }
private:
    static bool Nonempty(const int* start, const int* stop)
    {
        for (int n=0; n<3; ++n)
            if (start[n]<0 || start[n]>stop[n]) return false;
        return true;
    }
};
using Tiles3D = std::vector<Range3D>;

Tiles computeRectangularTilesNoDeps1D(int width, int blockWidth, int halfSteps);
Tiles computeDiamondTiles1D(int width, int blockWidth, int halfSteps);
std::vector<std::vector<Tiles3D>> computeRectangularTiles3D(int width[3], int blockWidth[3], int threads);
std::vector<std::vector<Tiles3D>> combineTilesTo3D(Tiles x, Tiles y, Tiles z, int halfSteps, int threads);

#endif

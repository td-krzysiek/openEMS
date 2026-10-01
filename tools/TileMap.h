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

#ifndef OPENEMS_TILE_MAP_H
#define OPENEMS_TILE_MAP_H
#include <vector>
#include <array>
#include <tuple>
#include <unordered_map>
#include <boost/functional/hash.hpp>

using TileKey = std::tuple<int, std::array<unsigned int, 3>, std::array<unsigned int, 3>>;

struct TileKeyHash {
	std::size_t operator()(const TileKey& key) const
	{
        std::size_t seed = 0;
        boost::hash_combine(seed, std::get<0>(key));
        for (unsigned int value : std::get<1>(key)) boost::hash_combine(seed, value);
        for (unsigned int value : std::get<2>(key)) boost::hash_combine(seed, value);
        return seed;
	}
};

inline TileKey GetTileKey(int order, const unsigned int start[3], const unsigned int stop[3])
{
	std::array<unsigned int, 3> startArray = {
		start[0], start[1], start[2]
	};
	std::array<unsigned int, 3> stopArray = {
		stop[0], stop[1], stop[2]
	};
	return std::make_tuple(order, startArray, stopArray);
}

using TileMap = std::unordered_map<TileKey, std::vector<int>, TileKeyHash>;

#endif

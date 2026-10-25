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

#include "tiling.h"
#include <stdexcept>

Tiles computeRectangularTilesNoDeps1D(int totalWidth, int blkWidth, int blkHalfTimesteps)
{
	if (totalWidth < 1 || blkWidth < 1 || blkHalfTimesteps != 2)
	{
		/*
		 * One timestep contains two half timesteps, one electric, one
		 * magnetic. Multi-timestep techniques (time skewing) are not
		 * supported.
		 */
		throw std::invalid_argument("Invalid tiling configuration");
	}

	std::vector<Block> blockList;

	/* at the beginning, divide the axis into several tiles */
	int numBlocks    = totalWidth / blkWidth;
	int numRemainder = totalWidth % blkWidth;

	/* For leftover cells, we allocate another block */
	if (numRemainder > 0) {
		numBlocks += 1;
	}

	blockList.resize(numBlocks);

	for (size_t i = 0; i < blockList.size(); i++) {
		Block block;

		for (int t = 0; t < 2; t++) {
			if (t % 2 == 0) {
				/* electric field range */
				Range startStop = Range(
					i * blkWidth,
					i * blkWidth + blkWidth - 1
				);
				if (startStop.second > totalWidth - 1) {
					startStop.second = totalWidth - 1;
				}
				block.push_back(startStop);
			}
			else {
				/* magnetic field range */
				Range startStop = Range(
					i * blkWidth,
					i * blkWidth + blkWidth - 1
				);
				if (startStop.second > totalWidth - 1) {
					startStop.second = totalWidth - 1;
				}
				block.push_back(startStop);
			}
		}

		blockList[i] = block;
	}

	/* TILES_RECTANGULAR has one phase. */
	Tiles tiles;
	tiles.type = TILES_RECTANGULAR;
	tiles.phases = 1;
	tiles.array.resize(tiles.phases);
	tiles.array[0] = blockList;

	return tiles;
}

std::vector<std::vector<Tiles3D>> computeRectangularTiles3D(
	int totalWidth[3],
	int blkWidth[3],
	int numThreads
)
{
	Tiles tilesX = computeRectangularTilesNoDeps1D(totalWidth[0], blkWidth[0], 2);
	Tiles tilesY = computeRectangularTilesNoDeps1D(totalWidth[1], blkWidth[1], 2);
	Tiles tilesZ = computeRectangularTilesNoDeps1D(totalWidth[2], blkWidth[2], 2);

	std::vector<std::vector<Tiles3D>> tilesPerPhasePerThread;
	tilesPerPhasePerThread.resize(numThreads);
	for (auto& tilesPerPhase : tilesPerPhasePerThread) {
		/* only one phase */
		tilesPerPhase.resize(1);
	}

	int assignedthread = 0;
	int blkHalfTimesteps = 2;

	for (int phaseX = 0; phaseX < tilesX.phases; phaseX++) {
		for (int phaseY = 0; phaseY < tilesY.phases; phaseY++) {
			for (int phaseZ = 0; phaseZ < tilesZ.phases; phaseZ++) {
				for (size_t x = 0; x < tilesX.array[phaseX].size(); x++) {
					for (size_t y = 0; y < tilesY.array[phaseY].size(); y++) {
						for (size_t z = 0; z < tilesZ.array[phaseZ].size(); z++) {
							for (int t = 0; t < blkHalfTimesteps; t += 2) {
								Range3D r;
								r.timestep = t / 2;

								r.voltageStart[0] = tilesX.array[phaseX][x][t].first;
								r.voltageStart[1] = tilesY.array[phaseY][y][t].first;
								r.voltageStart[2] = tilesZ.array[phaseZ][z][t].first;

								r.voltageStop[0]  = tilesX.array[phaseX][x][t].second;
								r.voltageStop[1]  = tilesY.array[phaseY][y][t].second;
								r.voltageStop[2]  = tilesZ.array[phaseZ][z][t].second;

								r.currentStart[0] = tilesX.array[phaseX][x][t + 1].first;
								r.currentStart[1] = tilesY.array[phaseY][y][t + 1].first;
								r.currentStart[2] = tilesZ.array[phaseZ][z][t + 1].first;

								r.currentStop[0]  = tilesX.array[phaseX][x][t + 1].second;
								r.currentStop[1]  = tilesY.array[phaseY][y][t + 1].second;
								r.currentStop[2]  = tilesZ.array[phaseZ][z][t + 1].second;

								bool omit = false;
								for (int n = 0; n < 3; n++) {
									if (r.voltageStart[n] == -1 && r.voltageStop[n] == -1) {
										omit = true;
									}
									if (r.currentStart[n] == -1 && r.currentStop[n] == -1) {
										omit = true;
									}
								}
								if (!omit) {
									tilesPerPhasePerThread[assignedthread][0].push_back(r);
									assignedthread = (assignedthread + 1) % numThreads;
								}
							}
						}
					}
				}
			}
		}
	}
	return tilesPerPhasePerThread;
}

Tiles computeDiamondTiles1D(int totalWidth, int blkWidth, int blkTimesteps)
{
	if (totalWidth < 1 || blkWidth < 1 || blkTimesteps < 2 || blkTimesteps % 2)
		throw std::invalid_argument("Invalid diamond tile dimensions");
	/*
	 * The width of a block is defined to be the shortest span
	 * of a mountain (top) or valley (bottom).
	 */
	int blkSpanMin = blkWidth;

	/*
	 * The longest span of a block is the mountain top or valley
	 * bottom.
	 */
	int blkSpanMax = blkWidth + blkTimesteps - 1;

	/* Calculate the total number of blocks at the last timestep. */
	int numBlocks    = totalWidth / (blkSpanMin + blkSpanMax) * 2;
	int numRemainder = totalWidth % (blkSpanMin + blkSpanMax);

	/* For leftover blocks, we add more mountains and valleys */
	for (int i = 0; numRemainder > 0; i++) {
		if (i % 2 == 0) {
			/* mountain */
			numBlocks++;
			numRemainder -= blkSpanMin;
		}
		else {
			/* valley */
			numBlocks++;
			numRemainder -= blkSpanMax;
		}
	}

	std::vector<Block> blockList;

	/*
	 * The rightmost block is possibly still incomplete (what a
	 * headache) and truncated on the edge. This happens when the
	 * last block is a full mountain, but its width decreases when
	 * we go back in time, creating empty space. Thus we always
	 * allocate one block more and pretend there's one more full
	 * block. We'll pretend the simulation space is large enough,
	 * and we remove out-of-bound blocks as the last step.
	 */
	blockList.resize(numBlocks + 1);

	/*
	 * Let's now consider how the axis should be splitted at the
	 * last timestep. It should be in the order of a mountain top
	 * (blkSpanMin), then a valley bottom (blkSpanMax), etc.
	 */
	int lastStop = -1;
	for (size_t i = 0; i < blockList.size(); i++) {
		Block block;
		block.resize(blkTimesteps);

		Range startStop;

		if (i % 2 == 0) {
			/* mountain */
			startStop = Range(
				lastStop + 1,
				lastStop + blkSpanMin
			);
		}
		else {
			/* valley */
			startStop = Range(
				lastStop + 1,
				lastStop + blkSpanMax
			);
		}
		lastStop = startStop.second;

		block[blkTimesteps - 1] = startStop;
		blockList[i] = block;
	}

	/* Finally we iterate timestamps in reverse, moving downwards. */
	for (int t = blkTimesteps - 2; t >= 0; t--) {
		for (size_t i = 0; i < blockList.size(); i++) {
			int prevStart = blockList[i][t + 1].first;
			int prevStop = blockList[i][t + 1].second;

			Range startStop;
			if (t % 2 != 0 && i % 2 == 0) {
				/*
				 * timestep is odd, block is a mountain,
				 * expand 1 unit to the left
				 */
				startStop = Range(
					prevStart - 1, prevStop
				);
			}
			else if (t % 2 == 0 && i % 2 == 0) {
				/*
				 * timestep is even, block is a mountain,
				 * expand 1 unit to the right
				 */
				startStop = Range(
					prevStart, prevStop + 1
				);
			}
			else if (t % 2 != 0 && i % 2 != 0) {
				/*
				 * timestep is odd, block is a valley,
				 * reduce 1 unit from right to the left
				 */
				startStop = Range(
					prevStart, prevStop - 1
				);
			}
			else if (t % 2 == 0 && i % 2 != 0) {
				/*
				 * timestep is even, block is a valley,
				 * reduce 1 unit from left to the right
				 */
				startStop = Range(
					prevStart + 1, prevStop
				);
			}

			/* truncate if off the left edge */
			if (startStop.first < 0) {
				startStop.first = 0;
			}
			blockList[i][t] = startStop;
		}
	}

	/*
	 * The headache of the last incomplete tile continues. We need
	 * to remove the parts that are outside the boundary.
	 */
	for (int t = 0; t < blkTimesteps; t++) {
		for (size_t i = 0; i < blockList.size(); i++) {
			Range r;
			r = blockList[i][t];
			if (r.first > totalWidth - 1) {
				r.first = -1;
				r.second = -1;
			}
			if (r.second > totalWidth - 1) {
				r.second = totalWidth - 1;
			}
			blockList[i][t] = r;
		}
	}

	/* TILES_DIAMOND has two phases. */
	Tiles tiles;
	tiles.type = TILES_DIAMOND;
	tiles.phases = 2;
	tiles.array.resize(tiles.phases);

	for (size_t i = 0; i < blockList.size(); i++) {
		if (i % 2 == 0) {
			/* all mountain blocks go to phase 1 */
			tiles.array[0].push_back(blockList[i]);
		}
		else {
			/* all valley blocks go to phase 2 */
			tiles.array[1].push_back(blockList[i]);
		}
	}

	return tiles;
}

std::vector<std::vector<Tiles3D>> combineTilesTo3D(
	Tiles tilesX, Tiles tilesY, Tiles tilesZ,
	int blkHalfTimesteps,
	int numThreads
)
{
	/*
	 * First, we need to find one dimension that is parallelizable
	 * (i.e. uses diamond tiling).
	 */
	int parallelAxis = -1;
	if (tilesZ.type == TILES_DIAMOND) {
		parallelAxis = 'Z';
	}
	else if (tilesY.type == TILES_DIAMOND) {
		parallelAxis = 'Y';
	}
	else if (tilesX.type == TILES_DIAMOND) {
		parallelAxis = 'X';
	}

	if (parallelAxis == -1 && numThreads != 1) {
		throw std::invalid_argument("Invalid tiling configuration");
	}

	int totalPhases = tilesX.phases * tilesY.phases * tilesZ.phases;
	std::vector<std::vector<Tiles3D>> tilesPerPhasePerThread;
	tilesPerPhasePerThread.resize(numThreads);

	for (auto& tilesPerPhase : tilesPerPhasePerThread) {
		tilesPerPhase.resize(totalPhases);
	}

	int phaseXYZ[3] = {0, 0, 0};  /* X, Y, Z */
	for (int phase = 0; phase < totalPhases; phase++) {
		int assignedThread = 0;
		int phaseX = phaseXYZ[0];
		int phaseY = phaseXYZ[1];
		int phaseZ = phaseXYZ[2];

		for (size_t x = 0; x < tilesX.array[phaseX].size(); x++) {
			for (size_t y = 0; y < tilesY.array[phaseY].size(); y++) {
				for (size_t z = 0; z < tilesZ.array[phaseZ].size(); z++) {
					for (int t = 0; t < blkHalfTimesteps; t += 2) {
						Range3D r;
						r.timestep = t / 2;

						r.voltageStart[0] = tilesX.array[phaseX][x][t].first;
						r.voltageStart[1] = tilesY.array[phaseY][y][t].first;
						r.voltageStart[2] = tilesZ.array[phaseZ][z][t].first;

						r.voltageStop[0]  = tilesX.array[phaseX][x][t].second;
						r.voltageStop[1]  = tilesY.array[phaseY][y][t].second;
						r.voltageStop[2]  = tilesZ.array[phaseZ][z][t].second;

						r.currentStart[0] = tilesX.array[phaseX][x][t + 1].first;
						r.currentStart[1] = tilesY.array[phaseY][y][t + 1].first;
						r.currentStart[2] = tilesZ.array[phaseZ][z][t + 1].first;

						r.currentStop[0]  = tilesX.array[phaseX][x][t + 1].second;
						r.currentStop[1]  = tilesY.array[phaseY][y][t + 1].second;
						r.currentStop[2]  = tilesZ.array[phaseZ][z][t + 1].second;

                        // A clipped tile can contain only E or only H in a
                        // half-step pair. Never drop the surviving update.
                        if (r.HasVoltage() || r.HasCurrent())
                            tilesPerPhasePerThread[assignedThread][phase].push_back(r);
					}
					if (parallelAxis == 'Z') {
						assignedThread = (assignedThread + 1) % numThreads;
					}
				}
				if (parallelAxis == 'Y') {
					assignedThread = (assignedThread + 1) % numThreads;
				}
			}
			if (parallelAxis == 'X') {
				assignedThread = (assignedThread + 1) % numThreads;
			}
		}

		/*
		 * Carry propagation for the (phaseX, phaseY, phaseZ) counter.
		 * This is equivelant to three nested for loops.
		 */
		phaseXYZ[2]++;

		if (phaseXYZ[2] > tilesZ.phases - 1) {
			phaseXYZ[2] = 0;
			phaseXYZ[1]++;
		}
		if (phaseXYZ[1] > tilesY.phases - 1) {
			phaseXYZ[1] = 0;
			phaseXYZ[0]++;
		}
		if (phaseXYZ[0] > tilesX.phases - 1) {
			/*
			 * No need to do anything if phaseX overflows, since the outer
			 * for loop will terminate when it reaches totalPhases.
			 */
		}
	}

	return tilesPerPhasePerThread;
}

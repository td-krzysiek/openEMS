/*
*	Copyright (C) 2023 Yifeng Li (tomli@tomli.me)
*	Copyright (C) 2010 Sebastian Held (sebastian.held@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "engine_tiling.h"
#include "extensions/engine_extension.h"
#include "tools/denormal.h"
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <memory>
#include <stdexcept>

namespace {
unsigned int tileSetting(const char* key, unsigned int fallback, unsigned int maximum)
{
    const char* value = std::getenv(key);
    if (!value) return fallback;
    char* end = NULL;
    unsigned long parsed = std::strtoul(value, &end, 10);
    if (!*value || *end || *value == '-' || parsed < 1 || parsed > maximum)
        throw std::invalid_argument(std::string(key) + " is outside the supported range");
    return static_cast<unsigned int>(parsed);
}
}

Engine_Tiling* Engine_Tiling::New(const Operator_Multithread* op, unsigned int numThreads)
{
    std::unique_ptr<Engine_Tiling> engine(new Engine_Tiling(op));
    engine->setNumThreads(numThreads);
    engine->Init();
    return engine.release();
}

Engine_Tiling::Engine_Tiling(const Operator_Multithread* op) : ENGINE_TILING_BASE(op),
    m_Op_MT(op), m_thread_group(NULL), m_startBarrier(NULL), m_stopBarrier(NULL),
    m_IterateBarrier(NULL), m_iterTS(0), m_numThreads(1), m_stopThreads(false),
    m_teamReady(false), m_startupFailed(false), m_tileX(10), m_tileY(10), m_tileSteps(5)
{
}

Engine_Tiling::~Engine_Tiling()
{
    Reset();
}

void Engine_Tiling::setNumThreads(unsigned int threads)
{
    m_numThreads = threads ? threads : std::max(1u, boost::thread::hardware_concurrency());
}

void Engine_Tiling::Init()
{
    ENGINE_TILING_BASE::Init();
    m_tileX = tileSetting("OPENEMS_TILE_X", 10, 4096);
    m_tileY = tileSetting("OPENEMS_TILE_Y", 10, 4096);
    m_tileSteps = tileSetting("OPENEMS_TILE_STEPS", 5, 64);
    changeNumThreads(m_numThreads);
}

void Engine_Tiling::Reset()
{
    if (m_thread_group)
    {
        // Public Reset is called between IterateTS batches. Release idle workers
        // normally; interrupting threads blocked at a barrier can deadlock.
        m_stopThreads = true;
        m_startBarrier->wait();
        m_thread_group->join_all();
        delete m_thread_group;
        m_thread_group = NULL;
    }
    delete m_startBarrier; m_startBarrier = NULL;
    delete m_stopBarrier; m_stopBarrier = NULL;
    delete m_IterateBarrier; m_IterateBarrier = NULL;
    ENGINE_TILING_BASE::Reset();
}

void Engine_Tiling::changeNumThreads(unsigned int threads)
{
    if (m_thread_group) throw std::logic_error("Tiling worker team is fixed during a run");
    int dimensions[3];
    for (unsigned int n=0; n<3; ++n)
    {
        if (numLines[n]<2 || numLines[n]>INT_MAX-128u)
            throw std::invalid_argument("Unsupported tiling grid size");
        dimensions[n] = static_cast<int>(numLines[n]);
    }
    int sizes[3] = {static_cast<int>(m_tileX), static_cast<int>(m_tileY), dimensions[2]};
    const int halfSteps = 2 * m_tileSteps;
    auto x = computeDiamondTiles1D(dimensions[0], sizes[0], halfSteps);
    auto y = computeDiamondTiles1D(dimensions[1], sizes[1], halfSteps);
    // SIMD packs the Z direction. Keep it whole and partition only XY.
    Tiles z;
    z.type = TILES_RECTANGULAR;
    z.phases = 1;
    z.array = {{Block(halfSteps, Range(0, dimensions[2] - 1))}};
    auto plan = combineTilesTo3D(x, y, z, halfSteps, threads);
    auto fallback = computeRectangularTiles3D(dimensions, sizes, threads);
    std::vector<Range3D> allTiles;
    for (const auto& team : plan)
        for (const auto& phase : team)
            allTiles.insert(allTiles.end(), phase.begin(), phase.end());
    for (const auto& team : fallback)
        for (const auto& phase : team)
            allTiles.insert(allTiles.end(), phase.begin(), phase.end());
    // Capability checks and all index-map allocations precede worker startup.
    InitializeTiling(allTiles);
    m_startBarrier = new boost::barrier(threads + 1);
    m_stopBarrier = new boost::barrier(threads + 1);
    m_IterateBarrier = new IterationBarrier(threads);
    m_stopThreads = false;
    m_thread_group = new boost::thread_group();
    m_teamReady = false;
    m_startupFailed = false;
    try
    {
        for (unsigned int n=0; n<threads; ++n)
            m_thread_group->create_thread(NS_Engine_Tiling::thread(this, plan[n], m_tileSteps, fallback[n], n));
    }
    catch (...)
    {
        { std::lock_guard<std::mutex> lock(m_startupMutex); m_startupFailed = true; }
        m_startupCondition.notify_all();
        m_thread_group->join_all();
        delete m_thread_group;
        m_thread_group = NULL;
        throw;
    }
    { std::lock_guard<std::mutex> lock(m_startupMutex); m_teamReady = true; }
    m_startupCondition.notify_all();
    std::cout << "Experimental tiling: " << m_tileX << "x" << m_tileY
              << "x(full Z), " << m_tileSteps << " steps, " << threads << " threads" << std::endl;
}

bool Engine_Tiling::IterateTS(unsigned int steps)
{
    if (steps > static_cast<unsigned int>(INT_MAX) - numTS)
        throw std::overflow_error("Tiling timestep exceeds signed index range");
    m_iterTS = steps;
    m_startBarrier->wait();
    m_stopBarrier->wait();
    return true;
}
void Engine_Tiling::DoPreVoltageUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in reverse order -> highest priority gets access to the voltages last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
	{
		m_Eng_exts.at(n)->DoPreVoltageUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}

}

void Engine_Tiling::DoPostVoltageUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->DoPostVoltageUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::Apply2Voltages(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the voltages first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->Apply2Voltages(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::DoPreCurrentUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in reverse order -> highest priority gets access to the currents last
	for (int n=m_Eng_exts.size()-1; n>=0; --n)
	{
		m_Eng_exts.at(n)->DoPreCurrentUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::DoPostCurrentUpdates(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->DoPostCurrentUpdates(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::Apply2Current(int timestep, int start[3], int stop[3])
{
	//execute extensions in normal order -> highest priority gets access to the currents first
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		m_Eng_exts.at(n)->Apply2Current(
			timestep,
			(unsigned int *) start,
			(unsigned int *) stop
		);
	}
}

void Engine_Tiling::InitializeTiling(const std::vector<Range3D>& tiles)
{
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		std::cerr << "Using ";
		std::cerr << m_Eng_exts.at(n)->GetExtensionName();
		std::cerr << std::endl;

		m_Eng_exts.at(n)->InitializeTiling(tiles);
	}
}

//
// *************************************************************************************************************************
//
namespace NS_Engine_Tiling
{

thread::thread(
	Engine_Tiling* ptr,
	std::vector<Tiles3D> tiles,
	unsigned int blkTimesteps,
	std::vector<Tiles3D> fallbackTiles,
	unsigned int threadID
)
{
	m_enginePtr = ptr;
	m_tiles = std::move(tiles);
	m_threadID = threadID;
	m_fallbackTiles = std::move(fallbackTiles);
	m_blkTimesteps = blkTimesteps;
}

void thread::operator()()
{
	Denormal::Disable();
    {
        std::unique_lock<std::mutex> lock(m_enginePtr->m_startupMutex);
        m_enginePtr->m_startupCondition.wait(lock, [this] {
            return m_enginePtr->m_teamReady || m_enginePtr->m_startupFailed;
        });
        if (m_enginePtr->m_startupFailed) return;
    }
	while (true)
	{
		m_enginePtr->m_startBarrier->wait();
		if (m_enginePtr->m_stopThreads) return;

		// how many timesteps are we calculating at the same time in
		// diamond tiling?
		unsigned int batchTimesteps = m_enginePtr->m_iterTS / m_blkTimesteps;
		unsigned int leftoverTimesteps = m_enginePtr->m_iterTS % m_blkTimesteps;

		int currentTimestep = m_enginePtr->numTS;

		for (unsigned int iter=0; iter < batchTimesteps; ++iter)
		{
			for (unsigned int stage = 0; stage < m_tiles.size(); stage++)
			{
				iterateTimesteps(currentTimestep, m_tiles[stage]);
				m_enginePtr->m_IterateBarrier->wait();
			}

			currentTimestep += m_blkTimesteps;
		}


		for (unsigned int iter = 0; iter < leftoverTimesteps; ++iter)
		{
			for (unsigned int stage = 0; stage < m_fallbackTiles.size(); stage++)
			{
				iterateUnskewedSingleTimestep(currentTimestep, m_fallbackTiles[stage]);
			}
			currentTimestep += 1;

		}

		if (m_threadID == 0) {
			 // only the first thread increments numTS
			m_enginePtr->numTS = currentTimestep;
		}
		m_enginePtr->m_stopBarrier->wait();
	}
}

void thread::updateVoltage(int timestep, Range3D& tile)
{
    if (!tile.HasVoltage()) return;
    unsigned int start[3], stop[3];
    for (int n=0; n<3; ++n) { start[n]=tile.voltageStart[n]; stop[n]=tile.voltageStop[n]; }
    m_enginePtr->DoPreVoltageUpdates(timestep, tile.voltageStart, tile.voltageStop);
    m_enginePtr->UpdateVoltagesTiled(start, stop);
    m_enginePtr->DoPostVoltageUpdates(timestep, tile.voltageStart, tile.voltageStop);
    m_enginePtr->Apply2Voltages(timestep, tile.voltageStart, tile.voltageStop);
}

void thread::updateCurrent(int timestep, Range3D& tile)
{
    if (!tile.HasCurrent()) return;
    unsigned int start[3], stop[3];
    for (int n=0; n<3; ++n)
    {
        start[n] = tile.currentStart[n];
        stop[n] = std::min(static_cast<unsigned int>(tile.currentStop[n]),
                          m_enginePtr->m_Op_MT->GetNumberOfLines(n) - 2);
    }
    m_enginePtr->DoPreCurrentUpdates(timestep, tile.currentStart, tile.currentStop);
    m_enginePtr->UpdateCurrentsTiled(start, stop);
    m_enginePtr->DoPostCurrentUpdates(timestep, tile.currentStart, tile.currentStop);
    m_enginePtr->Apply2Current(timestep, tile.currentStart, tile.currentStop);
}

void thread::iterateTimesteps(int timestep, std::vector<Range3D>& tiles)
{
    for (auto& tile : tiles)
    {
        updateVoltage(timestep + tile.timestep, tile);
        updateCurrent(timestep + tile.timestep, tile);
    }
}

void thread::iterateUnskewedSingleTimestep(int timestep, std::vector<Range3D>& tiles)
{
    for (auto& tile : tiles) updateVoltage(timestep, tile);
    m_enginePtr->m_IterateBarrier->wait();
    for (auto& tile : tiles) updateCurrent(timestep, tile);
    m_enginePtr->m_IterateBarrier->wait();
}

} // namespace

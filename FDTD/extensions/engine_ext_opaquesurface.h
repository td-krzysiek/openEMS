// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef ENGINE_EXT_OPAQUESURFACE_H
#define ENGINE_EXT_OPAQUESURFACE_H
#include "engine_extension.h"
#include <array>
#include <map>
#include <vector>
class Operator_Ext_OpaqueSurface;
class Engine_Ext_OpaqueSurface : public Engine_Extension
{
public:
    explicit Engine_Ext_OpaqueSurface(Operator_Ext_OpaqueSurface* op);
    virtual void InitializeTiling(const std::vector<Range3D>& tiles);
    virtual void Apply2Voltages();
    virtual void Apply2Voltages(int threadID);
    virtual void Apply2Voltages(int timestep, unsigned int start[3], unsigned int stop[3]);
private:
    template<class Eng> void Update(Eng* eng, const std::vector<size_t>& indices);
    Operator_Ext_OpaqueSurface* op;
    std::vector<double> voltage;
    std::vector<std::vector<double> > state, a, b;
    std::vector<size_t> all;
    std::map<std::array<unsigned int,6>, std::vector<size_t> > tile_map;
};
#endif

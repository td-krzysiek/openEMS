// SPDX-License-Identifier: GPL-3.0-or-later
#include "engine_ext_opaquesurface.h"
#include "operator_ext_opaquesurface.h"
#include "FDTD/engine_sse.h"
#include "engine_extension_dispatcher.h"

Engine_Ext_OpaqueSurface::Engine_Ext_OpaqueSurface(Operator_Ext_OpaqueSurface* p)
    : Engine_Extension(p), op(p)
{
    m_TilingSupported=true;
    // After PML, before sources. The operator rejects surface updates in PML.
    SetPriority(100);
    voltage.resize(op->edges.size(),0);
    state.resize(voltage.size()); a.resize(voltage.size()); b.resize(voltage.size());
    for (size_t i=0;i<voltage.size();++i) {
        all.push_back(i);
        const auto& e=op->edges[i]; const auto& f=op->fits[e.fit];
        state[i].resize(f.w.size(),0);
        for (size_t j=0;j<f.w.size();++j) {
            a[i].push_back((2*f.tau[j]-op->dt)/(2*f.tau[j]+op->dt));
            b[i].push_back(op->dt*f.w[j]*e.scale/(2*f.tau[j]+op->dt));
        }
    }
}

template<class Eng>
void Engine_Ext_OpaqueSurface::Update(Eng* eng, const std::vector<size_t>& indices)
{
    for (auto i : indices) {
        const auto& e=op->edges[i]; const unsigned int n=e.n;
        const unsigned int p=(n+1)%3, q=(n+2)%3;
        unsigned int lowp[3]={e.pos[0],e.pos[1],e.pos[2]};
        unsigned int lowq[3]={e.pos[0],e.pos[1],e.pos[2]};
        if (lowp[p]) --lowp[p];
        if (lowq[q]) --lowq[q];
        const double curl=eng->Eng::GetCurr(q,e.pos[0],e.pos[1],e.pos[2])
            -eng->Eng::GetCurr(q,lowp[0],lowp[1],lowp[2])
            -eng->Eng::GetCurr(p,e.pos[0],e.pos[1],e.pos[2])
            +eng->Eng::GetCurr(p,lowq[0],lowq[1],lowq[2]);
        double memory=0;
        for (size_t j=0;j<state[i].size();++j)
            memory+=(1+a[i][j])*state[i][j]/2;
        const double next=e.vv*voltage[i]+e.vi*(curl-memory);
        for (size_t j=0;j<state[i].size();++j)
            state[i][j]=a[i][j]*state[i][j]+b[i][j]*(next+voltage[i]);
        voltage[i]=next;
        eng->Eng::SetVolt(n,e.pos[0],e.pos[1],e.pos[2],next);
    }
}

void Engine_Ext_OpaqueSurface::Apply2Voltages() {ENG_DISPATCH_ARGS(Update,all);}
void Engine_Ext_OpaqueSurface::Apply2Voltages(int threadID)
{
    // One serial surface pass; main volume remains multithreaded.
    if (threadID==0) Apply2Voltages();
}

void Engine_Ext_OpaqueSurface::InitializeTiling(const std::vector<Range3D>& tiles)
{
    tile_map.clear();
    for (const auto& t: tiles) {
        if (!t.HasVoltage()) continue;
        std::array<unsigned int,6> key;
        for (int n=0;n<3;++n) {key[n]=t.voltageStart[n];key[n+3]=t.voltageStop[n];}
        auto entry=tile_map.emplace(key,std::vector<size_t>());
        if (!entry.second) continue;
        for (size_t i=0;i<op->edges.size();++i) {
            bool inside=true;
            for (int n=0;n<3;++n)
                inside=inside && op->edges[i].pos[n]>=key[n] && op->edges[i].pos[n]<=key[n+3];
            if (inside) entry.first->second.push_back(i);
        }
    }
}

void Engine_Ext_OpaqueSurface::Apply2Voltages(int, unsigned int start[3], unsigned int stop[3])
{
    const std::array<unsigned int,6> key={{start[0],start[1],start[2],stop[0],stop[1],stop[2]}};
    const auto& indices=tile_map.at(key);
    ENG_DISPATCH_ARGS(Update,indices);
}

# SPDX-License-Identifier: GPL-3.0-or-later
"""Shared fit coefficients must preserve the old per-edge FP32 bit patterns."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class OpaquePackingTest(unittest.TestCase):
    def test_shared_fits_independent_states_and_rounding(self):
        compiler = shutil.which(os.environ.get('CXX', 'c++'))
        if not compiler:
            self.skipTest('C++ compiler required for packing regression')
        root = Path(__file__).resolve().parents[2]
        source = r'''
#include "FDTD/gpu_opaquesurface.h"
#include <cassert>
#include <cstring>
struct Input {
    struct Fit { std::vector<double> w, tau; };
    struct Edge { unsigned n, pos[3], fit; double vv, vi, scale; };
    std::vector<Fit> fits;
    std::vector<Edge> edges;
    double dt;
};
bool bits(float x, float y) { return std::memcmp(&x,&y,sizeof(float))==0; }
int main() {
    Input op; op.dt=1.23456789e-14;
    for (unsigned order : {1u,62u,128u}) {
        Input::Fit f;
        for (unsigned j=0;j<order;++j) {
            f.w.push_back(.37+j*.013);
            f.tau.push_back(op.dt*(.27+j*.12));
        }
        op.fits.push_back(f);
    }
    for (unsigned i=0;i<12;++i) {
        Input::Edge e={i%3,{i%4,2,3},i%3,.75,.23,.1+i*.37};
        op.edges.push_back(e);
    }
    GPU_OpaqueSurfaceData d(op,9,8,7);
    assert(d.a.size()==191 && d.memory.size()==191);
    assert(d.b.size()==4*191 && d.offset.size()==13 && d.fit_offset.size()==12);
    bool rounding_sensitive=false;
    for (unsigned i=0;i<op.edges.size();++i) {
        const auto& e=op.edges[i]; const auto& f=op.fits[e.fit];
        assert(d.fit_offset[i]==d.fit_offset[i%3]);
        assert(d.offset[i+1]-d.offset[i]==f.w.size());
        assert(bits(d.vv[i],float(e.vv)) && bits(d.vi[i],float(e.vi)));
        const unsigned p=(e.n+1)%3,q=(e.n+2)%3;
        const uint64_t stride[3]={56,7,1}, cells=504;
        const uint64_t node=e.pos[0]*stride[0]+e.pos[1]*stride[1]+e.pos[2];
        const uint64_t expected[5]={e.n*cells+node,q*cells+node,
            q*cells+node-(e.pos[p]?stride[p]:0),p*cells+node,
            p*cells+node-(e.pos[q]?stride[q]:0)};
        for (unsigned k=0;k<5;++k) assert(d.index[5*i+k]==expected[k]);
        for (unsigned j=0;j<f.w.size();++j) {
            const double a=(2*f.tau[j]-op.dt)/(2*f.tau[j]+op.dt);
            assert(bits(d.a[d.fit_offset[i]+j],float(a)));
            assert(bits(d.memory[d.fit_offset[i]+j],float((1+a)/2)));
            rounding_sensitive |= !bits(float((1+a)/2),(1+float(a))/2);
            assert(bits(d.b[d.offset[i]+j],float(op.dt*f.w[j]*e.scale/(2*f.tau[j]+op.dt))));
        }
    }
    assert(rounding_sensitive); // catches recomputation from the rounded a
    assert(d.b[d.offset[0]]!=d.b[d.offset[3]]); // same fit, different edge scale
    // Distinct branch slices must not alias dynamic histories.
    std::vector<float> states(d.b.size(),0);
    for (unsigned i=0;i<op.edges.size();++i) states[d.offset[i]]=float(i+1);
    for (unsigned i=0;i<op.edges.size();++i) assert(states[d.offset[i]]==float(i+1));
    Input empty; empty.dt=op.dt;
    GPU_OpaqueSurfaceData z(empty,9,8,7);
    assert(z.a.empty() && z.b.empty() && z.fit_offset.empty() && z.offset.size()==1);
    bool rejected=false;
    op.edges[0].fit=99;
    try { GPU_OpaqueSurfaceData bad(op,9,8,7); } catch(const std::out_of_range&) {rejected=true;}
    assert(rejected); op.edges[0].fit=0;
    op.fits[0].tau.clear(); rejected=false;
    try { GPU_OpaqueSurfaceData bad(op,9,8,7); } catch(const std::runtime_error&) {rejected=true;}
    assert(rejected);
    for (double value : {std::numeric_limits<double>::infinity(),1e-300,1e300}) {
        rejected=false;
        try { GPU_OpaqueSurfaceData::Coefficient(value); } catch(const std::runtime_error&) {rejected=true;}
        assert(rejected);
    }
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            cpp = Path(tmp)/'test.cpp'
            cpp.write_text(source)
            binary = Path(tmp)/'test'
            subprocess.run([compiler, '-std=c++11', '-I', str(root), str(cpp),
                            '-o', str(binary)], check=True)
            subprocess.run([str(binary)], check=True)


if __name__ == '__main__':
    unittest.main()

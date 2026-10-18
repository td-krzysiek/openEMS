"""Standalone geometry regression: segmentation must not create metal faces."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class OpaqueUnionTest(unittest.TestCase):
    def test_segmentation_overlap_contacts_and_cavity(self):
        compiler=shutil.which(os.environ.get('CXX','c++'))
        if not compiler:
            self.skipTest('C++ compiler required for standalone union regression')
        root=Path(__file__).resolve().parents[2]
        source=r'''
#include "FDTD/extensions/opaque_surface_union.h"
#include <cassert>
void add(OpaqueSurfaceUnion& u,unsigned x0,unsigned y0,unsigned z0,
         unsigned x1,unsigned y1,unsigned z1,unsigned fit=0) {
    unsigned lo[3]={x0,y0,z0},hi[3]={x1,y1,z1};u.Add(lo,hi,fit);
}
int main() {
    OpaqueSurfaceUnion whole(8,8,8),split(8,8,8),overlap(8,8,8);
    add(whole,1,1,1,7,7,7);
    add(split,1,1,1,4,7,7);add(split,4,1,1,7,7,7);
    add(overlap,1,1,1,5,7,7);add(overlap,3,1,1,7,7,7);
    add(overlap,3,1,1,7,7,7); // duplicate must not add another sheet
    for (unsigned n=0;n<3;++n)
    for (unsigned x=0;x<=8;++x) for (unsigned y=0;y<=8;++y) for (unsigned z=0;z<=8;++z) {
        unsigned pos[3]={x,y,z};int a[2][2],b[2][2],c[2][2];
        whole.Quadrants(n,pos,a);split.Quadrants(n,pos,b);overlap.Quadrants(n,pos,c);
        for (int i=0;i<2;++i) for (int j=0;j<2;++j) assert(a[i][j]==b[i][j] && a[i][j]==c[i][j]);
        assert(whole.InteriorH(n,pos)==split.InteriorH(n,pos));
        assert(whole.InteriorH(n,pos)==overlap.InteriorH(n,pos));
    }
    unsigned seam[3]={4,3,3};assert(split.InteriorH(0,seam));
    OpaqueSurfaceUnion tee(8,8,8);
    add(tee,1,1,1,7,3,3);add(tee,3,3,1,5,7,3);
    unsigned junction[3]={3,3,2};int q[2][2];tee.Quadrants(2,junction,q);
    assert(q[0][0]>=0 && q[1][0]>=0 && q[1][1]>=0 && q[0][1]<0); // concave corner
    OpaqueSurfaceUnion cavity(8,8,8);
    add(cavity,1,1,1,3,7,7);add(cavity,5,1,1,7,7,7);
    add(cavity,3,1,1,5,3,7);add(cavity,3,5,1,5,7,7);
    assert(cavity.Get(4,4,4)<0 && cavity.Get(2,4,4)>=0); // preserve a hole
    bool rejected=false;try {add(whole,2,2,2,3,3,3,1);}catch(const std::runtime_error&) {rejected=true;}
    assert(rejected);
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            cpp=Path(tmp)/'test.cpp';cpp.write_text(source)
            binary=Path(tmp)/'test'
            subprocess.run([compiler,'-std=c++11','-I',str(root),str(cpp),'-o',str(binary)],check=True)
            subprocess.run([str(binary)],check=True)


if __name__=='__main__':unittest.main()

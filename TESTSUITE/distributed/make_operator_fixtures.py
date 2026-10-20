# SPDX-License-Identifier: GPL-3.0-or-later
"""Small static controls, not antenna simulations; no solver is launched."""
import argparse
from pathlib import Path
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    for graded in [False,True]:
        path=args.output/('graded.xml' if graded else 'uniform.xml')
        if path.exists():raise FileExistsError(path)
        csx=ContinuousStructure();grid=csx.GetGrid();grid.SetDeltaUnit(.001)
        for axis,n in zip('xyz',[35,19,13]):
            values=np.arange(n,dtype=float)*.1
            if graded:values+=.008*np.sin(np.arange(n)*.7)
            grid.AddLine(axis,values)
        material=csx.AddMaterial('dielectric',epsilon=3.556)
        material.AddBox([.8,.25,.2],[2.6,1.4,1.0],priority=1)
        metal=csx.AddMetal('PEC');metal.AddBox([.9,.9,.5],[2.4,1.3,.7],priority=2)
        fdtd=openEMS(NrTS=120,EndCriteria=0);fdtd.SetCSX(csx)
        fdtd.SetGaussExcite(24.125e9,2e9);fdtd.SetBoundaryCond(['PEC']*6)
        fdtd.Write2XML(str(path))

if __name__=='__main__':main()

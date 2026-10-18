# SPDX-License-Identifier: GPL-3.0-or-later
"""A passive opaque conductor must be independent of box segmentation.

CPU by default; --gpu additionally checks the device path. The 4096-step
window tests segmentation/engine equivalence, not energy convergence. This
isolated cube can retain a static excitation tail; the boundary is an RF limit.
"""
from pathlib import Path
import os
import sys
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path=Path(tempfile.gettempdir())/'Opaque_Surface_Union'
Sim_Path.mkdir(parents=True,exist_ok=True)
engines=['multithreaded']+(['gpu'] if '--gpu' in sys.argv else [])
reference=None
for engine in engines:
    for variant in ['whole','split','overlap','split_properties']:
        folder=Sim_Path/(engine+'-'+variant)
        fdtd=openEMS(NrTS=4096,EndCriteria=0)
        fdtd.SetGaussExcite(100e9,50e9);fdtd.SetBoundaryCond(['PML_4']*6)
        csx=ContinuousStructure();fdtd.SetCSX(csx)
        grid=csx.GetGrid();grid.SetDeltaUnit(1e-3)
        for axis in 'xyz':grid.AddLine(axis,np.linspace(0,2,21))
        metal=csx.AddMetal('opaque')
        metal.AddAttribute('ExperimentalOpaqueSurfaceV1','1;10,1e-12')
        spans=[(.8,1.2)] if variant=='whole' else [(.8,1),(1,1.2)] if variant.startswith('split') else [(.8,1.1),(.9,1.2),(.9,1.2)]
        for index,(lo,hi) in enumerate(spans):
            if variant=='split_properties' and index:
                metal=csx.AddMetal('opaque_second')
                metal.AddAttribute('ExperimentalOpaqueSurfaceV1','1;10,1e-12')
            metal.AddBox([lo,.8,.8],[hi,1.2,1.2])
        # A finite x slab avoids exact-plane roundoff in native CSXCAD output.
        csx.AddExcitation('source',exc_type=0,exc_val=[0,1,0]).AddBox([.55,.9,.9],[.65,1.1,1.1],priority=20)
        csx.AddProbe('voltage',p_type=0).AddBox([1.4,.8,.8],[1.4,.9,.8])
        sys.stdout.flush();saved=os.dup(1)
        with (Sim_Path/(engine+'-'+variant+'.log')).open('w') as log:
            try:
                os.dup2(log.fileno(),1)
                fdtd.Run(str(folder),cleanup=True,engine=engine,numThreads=2,exact_endcriteria=True,verbose=2)
            finally:
                sys.stdout.flush();os.dup2(saved,1);os.close(saved)
        text=(Sim_Path/(engine+'-'+variant+'.log')).read_text()
        assert 'Time for 4096 iterations' in text, 'Fixed regression window not completed'
        if engine=='gpu':
            assert 'all extensions run on the device' in text and 'host fallback' not in text
        data=np.loadtxt(folder/'voltage',comments='%',ndmin=2)
        assert np.isfinite(data).all() and np.max(abs(data[:,1]))>0, 'Inactive/nonfinite probe'
        if reference is None:reference=data
        assert data.shape==reference.shape, 'Different stopping/sampling schedule'
        error=np.max(abs(data[:,1]-reference[:,1]))/np.max(abs(reference[:,1]))
        tolerance=5e-4 if engine=='gpu' else 1e-6
        # Geometry permutations share the exact EC coefficients. The GPU
        # allowance covers FP32/FMA ordering relative to CPU boundary state.
        assert error<tolerance, 'FAIL: relative peak error {} >= {}'.format(error,tolerance)
        print(engine,variant,'relative peak error',error,'tolerance',tolerance)
print('PASS')

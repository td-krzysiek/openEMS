"""Native compact device reconstruction/probes/FD dumps against the double CPU reference.

Default: explicitly selected gpu-reference, usable on CPU CI. Set
OPENEMS_SUBGRID_DEVICE=metal or cuda for a separately checked physical build.
No backend/device availability failure is skipped in the requested device run.
"""
import json
import sys
import os
import tempfile
from pathlib import Path

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from _subgrid_device_support import run_native_device

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_Device')
os.makedirs(Sim_Path, exist_ok=True)
PHYSICAL = os.environ.get('OPENEMS_SUBGRID_DEVICE', '')
assert PHYSICAL in ('', 'metal', 'cuda'), 'unsupported qualification selector'
GATE = 2e-4


def fixture(kind='local', steps=2048, source_x=-10):
    fdtd = openEMS(NrTS=steps, EndCriteria=0, TimeStepFactor=0.8)
    if kind == 'multiple':
        fdtd.SetCartesianSubgrids([
            {'box': [-8, -4, -4, -4, 4, 4], 'ratio': [2, 2, 4]},
            {'box': [4, -4, -4, 8, 4, 4], 'ratio': [4, 1, 2]},
        ])
    else:
        fdtd.SetCartesianSubgrid([-4]*3+[4]*3)
        fdtd.SetCartesianSubgridRatio([2, 1, 4] if kind == 'graded' else [2]*3)
    fdtd.SetGaussExcite(13e9, 9e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis in 'xyz':
        lines = np.arange(-16, 17, 2, dtype=float)
        if kind == 'graded':
            # Preserve exact region boundaries and source/probe coordinates.
            lines[lines == -2] = -2.3
            lines[lines == 6] = 6.3
        grid.SetLines(axis, lines)
    csx.AddExcitation('source', exc_type=0, exc_val=[0, 0, 1]).AddBox([source_x, 0, 0], [source_x, 0, 2])
    csx.AddProbe('electric', p_type=2).AddPoint([2, 2, 2])
    csx.AddProbe('magnetic', p_type=3).AddPoint([2, 2, 2])
    for name, quantity in [('ef', 10), ('hf', 11)]:
        csx.AddDump(name, dump_type=quantity, file_type=1, frequency=[10e9,13e9,16e9]).AddBox([-4,-4,0], [4,4,0])
    if kind == 'material-pec':
        csx.AddMaterial('dielectric', epsilon=3.556).AddBox([0,-16,-16], [16,16,16])
        csx.AddMetal('pec-sheet').AddBox([-8,-4,-4], [-8,4,4])
    return fdtd, csx


def run(kind, engine, **options):
    fdtd, csx = fixture(kind)
    previous = os.getcwd()
    path = os.path.join(Sim_Path, kind+'-'+engine)
    try:
        run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1, exact_endcriteria=True, **options)
    finally:
        os.chdir(previous)
    return Path(path)


def compare(reference, candidate):
    errors = {}
    for name in ('electric','magnetic'):
        a, b = [np.loadtxt(path/name, comments='%') for path in (reference,candidate)]
        assert a.shape == b.shape and np.isfinite(b).all()
        np.testing.assert_allclose(a[:,0],b[:,0],rtol=1e-12,atol=0)
        errors[name] = float(np.max(abs(a[:,1:]-b[:,1:]))/np.max(abs(a[:,1:])))
    for name in ('ef.h5','hf.h5'):
        with h5py.File(reference/name) as a, h5py.File(candidate/name) as b:
            keys = []
            a['FieldData'].visititems(lambda key,obj: keys.append(key) if isinstance(obj,h5py.Dataset) else None)
            assert keys
            scale = max(float(np.max(abs(a['FieldData/'+key][...]))) for key in keys)
            delta = 0
            for key in keys:
                left, right = a['FieldData/'+key][...], b['FieldData/'+key][...]
                assert left.shape == right.shape and np.isfinite(right).all()
                delta = max(delta,float(np.max(abs(left-right))))
            errors[name] = delta/scale
    assert all(error < GATE for error in errors.values()), (candidate,errors)
    print(candidate.name, errors)
    return errors


results = {}
for kind in ('local','graded','multiple','material-pec'):
    baseline = run(kind,'basic')
    reference = run(kind,'gpu-reference')
    results[kind+'-reference'] = compare(baseline,reference)
    if PHYSICAL:
        device = run(kind,'gpu',verbose=1)
        results[kind+'-'+PHYSICAL] = compare(baseline,device)
        results[kind+'-'+PHYSICAL+'-v-reference'] = compare(reference,device)

# Native engine construction errors return setup code4 through the current Python API.
# Capture native stderr so an unrelated setup failure cannot satisfy this control.
for kind in ('mur',):
    fdtd, csx = fixture(steps=8, source_x=0)
    if kind == 'resistor':
        fdtd.AddLumpedPort(1,50,[-10,0,0],[-10,0,2],'z',1)
    else:
        fdtd.SetBoundaryCond(['MUR','MUR','PEC','PEC','PEC','PEC'])
    previous = os.getcwd()
    saved = os.dup(2)
    guard_log = Path(Sim_Path,'reject-'+kind+'.stderr')
    try:
        with guard_log.open('w') as log:
            sys.stderr.flush()
            os.dup2(log.fileno(),2)
            result = fdtd.Run(os.path.join(Sim_Path,'reject-'+kind),
                              engine='gpu' if PHYSICAL else 'gpu-reference',numThreads=1,
                              gpu_device=0 if PHYSICAL == 'cuda' else -1)
            assert result == 4, (kind,result)
    finally:
        os.dup2(saved,2)
        os.close(saved)
        os.chdir(previous)
    reason = guard_log.read_text()
    expected = 'device resistor coupling' if kind == 'resistor' else 'device absorber/extension'
    assert 'CartesianSubgrid' in reason and expected in reason, reason
    assert not Path(Sim_Path,'reject-'+kind,'electric').exists(), 'unsupported case stepped'
    print('correctly rejected',kind,reason.strip())
Path(Sim_Path,'metrics.json').write_text(json.dumps({'gate':GATE,'physical':PHYSICAL,'errors':results},indent=2)+'\n')
print('PASS: native compact device E/H and FD reconstruction; scoped unsupported-physics guards')

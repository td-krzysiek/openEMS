"""Exact child-coordinate API/XML propagation and complete-batch E/H controls."""
from pathlib import Path
import json
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

ROOT = Path(tempfile.gettempdir()) / 'Cartesian_Subgrid_ChildLines'
ROOT.mkdir(exist_ok=False)

def fixture():
    base = np.arange(-8., 9., 2.)
    children = []
    for n in range(3):
        values = []
        for lo, hi in zip(base[:-1], base[1:]):
            values.extend([lo, lo + (.25 + .125*n)*(hi-lo)])
        children.append(values + [base[-1]])
    fdtd = openEMS(NrTS=2048, EndCriteria=0, TimeStepFactor=.8,
                  CartesianSubgridBox=[-4,-4,-4,4,4,4],
                  CartesianSubgridChildLines=children)
    fdtd.SetGaussExcite(13e9, 9e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(.001)
    for axis in 'xyz':
        mesh.SetLines(axis, base)
    csx.AddExcitation('source', exc_type=0, exc_val=[0,0,1]).AddBox([-6,0,-2], [-6,0,2])
    csx.AddProbe('received', p_type=2).AddPoint([2,2,2])
    csx.AddProbe('magnetic', p_type=3).AddPoint([2,2,2])
    return fdtd, children

def run(fdtd, name, engine='basic', setup_only=False):
    previous = os.getcwd()
    try:
        fdtd.Run(str(ROOT/name), engine=engine, setup_only=setup_only,
                 numThreads=1, exact_endcriteria=True)
    finally:
        os.chdir(previous)

first, axes = fixture()
run(first, 'api-basic')
first.Write2XML(str(ROOT/'roundtrip.xml'))
restored = openEMS()
assert restored.ReadFromXML(str(ROOT/'roundtrip.xml'))
run(restored, 'xml-basic')
device, _ = fixture()
run(device, 'api-reference', engine='gpu-reference')
metrics = []
for quantity in ['received', 'magnetic']:
    expected = np.loadtxt(ROOT/'api-basic'/quantity, comments='%')
    assert len(expected) > 10 and np.isfinite(expected).all()
    scale = np.max(np.abs(expected[:,1:]))
    assert scale > 1e-7, (quantity, scale)
    for case in ['xml-basic', 'api-reference']:
        actual = np.loadtxt(ROOT/case/quantity, comments='%')
        assert np.array_equal(actual[:,0], expected[:,0])
        error = float(np.max(np.abs(actual[:,1:] - expected[:,1:]))/scale)
        assert error < (1e-14 if case == 'xml-basic' else 2e-4), (case, quantity, error)
        metrics.append(dict(case=case, quantity=quantity, relativeError=error))

for name, change in [('bad-anchor', lambda a: a[0].__setitem__(4, a[0][4]+.1)),
                     ('bad-length', lambda a: a[1].pop())]:
    invalid, children = fixture()
    change(children)
    invalid.SetCartesianSubgridChildLines(children)
    try:
        run(invalid, name, setup_only=True)
    except (ValueError, RuntimeError) as error:
        assert 'CartesianSubgrid' in str(error), str(error)
    else:
        raise AssertionError(name + ' was accepted')

result = dict(status='PASS_NATIVE_CHILD_COORDINATE_API_XML_EH', validCases=3,
              negativeMeshCases=2, observations=metrics,
              limits=['Public vacuum fixture only; no PCB RF, Metal/CUDA or scaling qualification.'])
(ROOT/'metrics.json').write_text(json.dumps(result, indent=2)+'\n')
print(result['status'])

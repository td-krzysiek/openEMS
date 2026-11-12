"""Native planar-grid setup, source/probe integration and fail-closed guards.

The C++ CTest covers numerical convergence and late-time stability. Here a
small empty PEC cavity exercises the actual Python -> solver -> probe path.
"""
import os
import tempfile

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


def fixture(**options):
    fdtd = openEMS(NrTS=1024, EndCriteria=0, TimeStepFactor=0.8, CartesianMultiGridX=0, **options)
    fdtd.SetGaussExcite(8e9, 6e9)
    fdtd.SetBoundaryCond(['PEC'] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    grid.SetLines('x', np.arange(-16, 17))
    grid.SetLines('y', np.arange(-16, 17))
    grid.SetLines('z', np.arange(-4, 5))
    csx.AddExcitation('source', exc_type=0, exc_val=[0, 0, 1]).AddBox([-8, 0, -1], [-8, 0, 1])
    csx.AddProbe('received', p_type=2).AddPoint([8, 0, 0])
    return fdtd, csx


def run(fdtd, name, **options):
    previous = os.getcwd()
    try:
        return fdtd.Run(os.path.join(Sim_Path, name), numThreads=1, exact_endcriteria=True, **options)
    finally:
        os.chdir(previous)


def rejected(name, change, **options):
    fdtd, csx = fixture()
    change(fdtd, csx)
    try:
        run(fdtd, name, setup_only=True, engine=options.get('engine', 'basic'))
    except (ValueError, RuntimeError) as error:
        assert 'CartesianMultiGrid' in str(error), str(error)
        print(name + ': correctly rejected: ' + str(error))
        return
    assert False, name + ': unsupported setup was silently accepted'


Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_MultiGrid')
os.makedirs(Sim_Path, exist_ok=True)
fdtd, csx = fixture()
run(fdtd, 'valid', engine='basic')
received = np.loadtxt(os.path.join(Sim_Path, 'valid', 'received'), comments='%')
assert received.shape[0] > 10 and np.isfinite(received).all(), 'invalid native probe recording'
assert np.max(np.abs(received[:, 1:])) > 1e-5, 'no wave crossed the interface'

rejected('even-y', lambda f, c: c.GetGrid().SetLines('y', np.arange(-16, 16)))
rejected('nonuniform-y', lambda f, c: c.GetGrid().SetLines('y', np.arange(33) ** 1.1))
rejected('unanchored', lambda f, c: f.SetCartesianMultiGrid(0.5))
rejected('pml', lambda f, c: f.SetBoundaryCond(['PML_2'] * 6))
rejected('gpu', lambda f, c: None, engine='gpu')
rejected('cylindrical', lambda f, c: f.SetCylinderCoords())
rejected('crossing-metal', lambda f, c: c.AddMetal('crossing').AddBox([-1, -1, -1], [1, 1, 1]))
rejected('crossing-material', lambda f, c: c.AddMaterial('crossing', epsilon=2).AddBox([-1, -1, -1], [1, 1, 1]))
rejected('crossing-probe', lambda f, c: c.AddProbe('bad', p_type=3).AddPoint([0, 0, 0]))
rejected('tfsf', lambda f, c: c.AddExcitation('plane', exc_type=10, exc_val=[0, 0, 1]).AddBox([-12, -2, -2], [-10, 2, 2]))
rejected('lumped-rc', lambda f, c: c.AddLumpedElement('rc', ny='z', R=50, C=1e-12).AddBox([-10, 0, -1], [-10, 0, 1]))
rejected('transformed', lambda f, c: c.AddMetal('translated').AddBox([-10, -1, -1], [-9, 1, 1]).GetTransform().AddTransform('Translate', [10, 0, 0]))
print('PASS: Cartesian multi-grid source/probe integration and unsupported-case rejection')

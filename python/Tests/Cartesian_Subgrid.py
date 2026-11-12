"""Native compact XYZ run, E/H processing and fail-closed scope checks.

Numerical stability, native Yee equivalence and convergence are CTest tests.
Here both a probe and a field dump read through the compact engine adapter.
"""
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


def fixture(ratio=None):
    fdtd = openEMS(NrTS=2048, EndCriteria=0, TimeStepFactor=0.8,
                  CartesianSubgridBox=[-4, -4, -4, 4, 4, 4])
    if ratio is not None:
        fdtd.SetCartesianSubgridRatio(ratio)
    fdtd.SetGaussExcite(13e9, 9e9)
    fdtd.SetBoundaryCond(['PEC'] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis in 'xyz':
        grid.SetLines(axis, np.arange(-8, 9, 2))
    csx.AddExcitation('source', exc_type=0, exc_val=[0, 0, 1]).AddBox([-6, 0, -2], [-6, 0, 2])
    csx.AddProbe('received', p_type=2).AddPoint([2, 2, 2])
    csx.AddProbe('magnetic', p_type=3).AddPoint([2, 2, 2])
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
        assert 'CartesianSubgrid' in str(error), str(error)
        print(name + ': correctly rejected: ' + str(error))
        return
    assert False, name + ': unsupported setup silently accepted'


Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid')
os.makedirs(Sim_Path, exist_ok=True)
fdtd, csx = fixture()
csx.AddDump('ef', dump_type=10, file_type=1, frequency=[13e9]).AddBox([-4, -4, 0], [4, 4, 0])
csx.AddDump('hf', dump_type=11, file_type=1, frequency=[13e9]).AddBox([-4, -4, 0], [4, 4, 0])
run(fdtd, 'valid', engine='basic')
for name in ('received', 'magnetic'):
    data = np.loadtxt(os.path.join(Sim_Path, 'valid', name), comments='%')
    assert data.shape[0] > 10 and np.isfinite(data).all(), 'invalid native probe data'
    assert np.max(np.abs(data[:, 1:])) > 1e-7, 'no propagated field at probe'
for name in ('ef.h5', 'hf.h5'):
    assert os.path.getsize(os.path.join(Sim_Path, 'valid', name)) > 100, 'field dump missing'

# Exercise actual adapter propagation and E/H processing with unequal ratios.
fdtd, csx = fixture([1, 2, 4])
run(fdtd, 'unequal', engine='basic')
for name in ('received', 'magnetic'):
    data = np.loadtxt(os.path.join(Sim_Path, 'unequal', name), comments='%')
    assert np.isfinite(data).all() and np.max(np.abs(data[:, 1:])) > 1e-7

fdtd, csx = fixture()
fdtd.SetCartesianSubgrids([
    {'box': [-6, -4, -4, -2, 4, 4], 'ratio': [2, 2, 4]},
    {'box': [2, -4, -4, 6, 4, 4], 'ratio': [4, 1, 2]},
])
csx.AddDump('multi_ef', dump_type=10, file_type=1, frequency=[13e9]).AddBox([-6, -4, 0], [6, 4, 0])
csx.AddDump('multi_hf', dump_type=11, file_type=1, frequency=[13e9]).AddBox([-6, -4, 0], [6, 4, 0])
run(fdtd, 'multiple', engine='basic')
for name in ('received', 'magnetic'):
    data = np.loadtxt(os.path.join(Sim_Path, 'multiple', name), comments='%')
    assert np.isfinite(data).all() and np.max(np.abs(data[:, 1:])) > 1e-7
for name in ('multi_ef.h5', 'multi_hf.h5'):
    assert os.path.getsize(os.path.join(Sim_Path, 'multiple', name)) > 100
rejected('multi-unaligned', lambda f,c: f.SetCartesianSubgrids([
    {'box': [-6, -4, -4, -2.1, 4, 4]}, {'box': [2, -4, -4, 6, 4, 4]}]))
rejected('multi-boundary', lambda f,c: f.SetCartesianSubgrids([
    {'box': [-8, -4, -4, -2, 4, 4]}, {'box': [2, -4, -4, 6, 4, 4]}]))
rejected('multi-gpu', lambda f,c: f.SetCartesianSubgrids([
    {'box': [-6, -4, -4, -2, 4, 4]}, {'box': [2, -4, -4, 6, 4, 4]}]), engine='gpu')

rejected('pml', lambda f, c: f.SetBoundaryCond(['PML_2'] * 6))
rejected('sse', lambda f, c: None, engine='sse')
rejected('gpu', lambda f, c: None, engine='gpu')
rejected('nonuniform', lambda f, c: c.GetGrid().SetLines('z', np.arange(9) ** 1.1))
rejected('unaligned', lambda f, c: f.SetCartesianSubgrid([-3.9, -4, -4, 4, 4, 4]))
rejected('boundary-box', lambda f, c: f.SetCartesianSubgrid([-8, -4, -4, 4, 4, 4]))
rejected('material', lambda f, c: c.AddMaterial('slab', epsilon=2).AddBox([-1]*3, [1]*3))
rejected('metal', lambda f, c: c.AddMetal('patch').AddBox([-1]*3, [1]*3))
rejected('port', lambda f, c: f.AddLumpedPort(1, 50, [-6, 0, 0], [-6, 0, 2], 'z', 1))
rejected('bad-step', lambda f, c: f.SetTimeStep(1e-8))
rejected('d-field', lambda f, c: c.AddDump('d', dump_type=4).AddBox([-2]*3, [2]*3))
rejected('hard-source', lambda f, c: c.AddExcitation('hard', exc_type=1, exc_val=[0, 0, 1]).AddBox([-6, 0, 0], [-6, 0, 2]))
rejected('empty-source', lambda f, c: c.AddExcitation('missed', exc_type=0, exc_val=[0, 0, 1]).AddBox([-6.7, 0.3, 0], [-6.7, 0.3, 2]))
print('PASS: compact XYZ native run, E/H probes/dumps and unsupported-case rejection')

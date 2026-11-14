"""Lossless cell-aligned materials: native priority, units and strict guards.

Analytic layered cavity, heterogeneous invariant and native-Yee comparison
are CTest controls. These tests exercise the scripting/native mapping path.
"""
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_Materials')
os.makedirs(Sim_Path, exist_ok=True)


def fixture(offset=(0, 0, 0), drawing_scale=1, reverse=False):
    offset = np.asarray(offset)

    def point(value):
        return ((np.asarray(value) + offset) * drawing_scale).tolist()

    fdtd = openEMS(NrTS=1024, EndCriteria=0, TimeStepFactor=0.8,
                  CartesianSubgridBox=point([-4]*3) + point([4]*3))
    fdtd.SetGaussExcite(13e9, 9e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3/drawing_scale)
    for n, axis in enumerate('xyz'):
        grid.SetLines(axis, ((np.arange(-8, 9, 2)+offset[n])*drawing_scale).tolist())
    csx.AddExcitation('source', exc_type=0, exc_val=[0, 0, 1]).AddBox(
        point([-6, 0, -2]), point([-6, 0, 2]))
    csx.AddProbe('received', p_type=2).AddPoint(point([2, 2, 2]))
    csx.AddProbe('magnetic', p_type=3).AddPoint(point([2, 2, 2]))
    definitions = [('slab', 3.556, 1, 2, [0, -8, -8], [8, 8, 8]),
                   ('priority', 6, 1.7, 9, [0, -2, -2], [2, 2, 2])]
    for name, eps, mue, priority, start, stop in (reversed(definitions) if reverse else definitions):
        csx.AddMaterial(name, epsilon=eps, mue=mue).AddBox(point(start), point(stop), priority=priority)
    return fdtd, csx


def run(fdtd, label, setup_only=False):
    previous = os.getcwd()
    try:
        fdtd.Run(os.path.join(Sim_Path, label), engine='basic', numThreads=1,
                 exact_endcriteria=True, setup_only=setup_only)
    finally:
        os.chdir(previous)


for label, options in [('reference', {}), ('reverse', {'reverse': True}),
                       ('translated', {'offset': [31, -23, 17]}),
                       ('micrometre-units', {'drawing_scale': 1000})]:
    fdtd, csx = fixture(**options)
    run(fdtd, label)
    for name in ['received', 'magnetic']:
        data = np.loadtxt(os.path.join(Sim_Path, label, name), comments='%')
        assert np.isfinite(data).all() and np.max(np.abs(data[:, 1:])) > 1e-7
        if label != 'reference':
            reference = np.loadtxt(os.path.join(Sim_Path, 'reference', name), comments='%')
            assert data.shape == reference.shape
            np.testing.assert_allclose(data[:, 0], reference[:, 0], rtol=1e-11, atol=0)
            error = np.max(np.abs(data[:, 1:] - reference[:, 1:])) / np.max(np.abs(reference[:, 1:]))
            print(label, name, 'peak-relative error', error)
            assert error < 1e-8, (label, name, error)

# A material may be aligned to local fine cells without being on base lines.
fdtd, csx = fixture()
csx.AddMaterial('local', epsilon=2.25).AddBox([-3, -3, -3], [-1, -1, -1], priority=11)
run(fdtd, 'local-fine-alignment')


def reject(label, change):
    fdtd, csx = fixture()
    change(csx)
    try:
        run(fdtd, label, setup_only=True)
    except (ValueError, RuntimeError) as error:
        assert 'CartesianSubgrid' in str(error), str(error)
        print(label, 'correctly rejected:', error)
        return
    assert False, label + ': unsupported material silently accepted'


reject('unaligned', lambda c: c.AddMaterial('bad', epsilon=2).AddBox([-3.1]*3, [-1]*3))
reject('almost-aligned', lambda c: c.AddMaterial('bad', epsilon=2).AddBox([-3+1e-9]*3, [-1]*3))
reject('zero-volume', lambda c: c.AddMaterial('bad', epsilon=2).AddBox([-2]*3, [-2, 2, 2]))
reject('negative', lambda c: c.AddMaterial('bad', epsilon=-1).AddBox([-2]*3, [2]*3))
reject('lossy-electric', lambda c: c.AddMaterial('bad', epsilon=2, kappa=.1).AddBox([-2]*3, [2]*3))
reject('lossy-magnetic', lambda c: c.AddMaterial('bad', epsilon=2, sigma=.1).AddBox([-2]*3, [2]*3))
reject('non-box', lambda c: c.AddMaterial('bad', epsilon=2).AddSphere([0, 0, 0], 2))
reject('non-cartesian', lambda c: c.AddMaterial('bad', epsilon=2).AddBox([-2]*3, [2]*3).SetCoordinateSystem(1))
reject('transformed', lambda c: c.AddMaterial('bad', epsilon=2).AddBox([-2]*3, [2]*3).AddTransform('Translate', [1, 0, 0]))


def anisotropic(csx):
    material = csx.AddMaterial('bad', epsilon=[2, 3, 4])
    material.SetIsotropy(False)
    material.AddBox([-2]*3, [2]*3)


def weighted(csx, function):
    material = csx.AddMaterial('bad', epsilon=2)
    material.SetMaterialWeight(epsilon=function)
    material.AddBox([-2]*3, [2]*3)


reject('anisotropic', anisotropic)
reject('weight-spatial', lambda c: weighted(c, '1+x*x'))
reject('weight-constant', lambda c: weighted(c, '2'))

# An otherwise distinct base mesh may have unrepresentable child centres at
# an enormous drawing origin. Do not classify those leaves at rounded points.
fdtd, csx = fixture(offset=[1e16]*3)
try:
    run(fdtd, 'unrepresentable-origin', setup_only=True)
except (ValueError, RuntimeError) as error:
    assert 'CartesianSubgrid' in str(error), str(error)
else:
    assert False, 'unrepresentable native material sampling was silently accepted'
print('PASS: native subgrid material priority, coordinates, units, local alignment and guards')

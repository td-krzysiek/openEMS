# SPDX-License-Identifier: GPL-3.0-or-later
"""PMC symmetry planes against independently reflected closed PEC domains.

The compact reference places PMC at the physical outer base face. Ordinary
Yee PMC lives on staggered magnetic planes, so directly reusing its boundary
line coordinates would change the physical cavity. Instead use ordinary Yee
on the full PEC domain with mirror-related sources as the identity control.
Refined half domains compare with the matched full two-region compact grid.
"""
import json
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_Boundaries')
os.makedirs(Sim_Path, exist_ok=True)
metrics = {}


def run(normal, positive, half, refined, graded=False):
    label = '{}{}-{}-{}'.format('xyz'[normal], '+' if positive else '-',
                              'half' if half else 'full', ('local-graded' if graded else 'local') if refined else 'identity')
    kwargs = dict(NrTS=3072, EndCriteria=0, TimeStep=5e-13, OverSampling=10)
    fdtd = openEMS(**kwargs)
    fdtd.SetGaussExcite(13e9, 9e9)
    bc = ['PEC']*6
    if half:
        bc[2*normal + (0 if positive else 1)] = 'PMC'
    fdtd.SetBoundaryCond(bc)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for n, axis in enumerate('xyz'):
        lines = np.array([-8,-6.3,-4,-2,0,2,4,6.3,8]) if graded else np.arange(-8,9,2)
        if n == normal and half:
            lines = lines[lines >= 0] if positive else lines[lines <= 0]
        grid.SetLines(axis, lines)
    boxes = []
    for side in ([1 if positive else -1] if half else [-1, 1]):
        lo = [-4.]*3
        hi = [4.]*3
        outer = 6.3 if graded else 6.
        lo[normal], hi[normal] = ((2., outer) if side == 1 else (-outer, -2.))
        boxes.append([lo, hi])
    if half or refined:
        ratio = [2, 1, 2] if refined else [1, 1, 1]
        for lo, hi in boxes:
            fdtd.AddCartesianSubgrid(box=lo+hi, ratio=ratio)
    source_axis = (normal+1) % 3
    vector = [0., 0., 0.]
    vector[source_axis] = 1.
    for side in ([1 if positive else -1] if half else [-1, 1]):
        lo = [0.]*3
        hi = [0.]*3
        lo[normal] = hi[normal] = side*4.
        lo[source_axis], hi[source_axis] = -2., 2.
        csx.AddExcitation('source'+str(side), exc_type=0, exc_val=vector).AddBox(lo, hi)
    # Three matched interior points, plus tangential E on the symmetry plane.
    for i in range(4):
        point = [2., 2., 2.]
        point[normal] = ((1 if positive else -1) * [0., 2., 4., 6.][i])
        csx.AddProbe('e'+str(i), p_type=2).AddPoint(point)
        csx.AddProbe('h'+str(i), p_type=3).AddPoint(point)
    path = os.path.join(Sim_Path, label)
    os.makedirs(path, exist_ok=True)
    fdtd.Write2XML(os.path.join(path, 'input.xml'))
    old = os.getcwd()
    try:
        fdtd.Run(path, engine='basic', numThreads=1, exact_endcriteria=True)
    finally:
        os.chdir(old)
    return {name:np.loadtxt(os.path.join(path, name), comments='%')
            for name in ['e'+str(i) for i in range(4)]+['h'+str(i) for i in range(4)]}


for refined, graded in [(False, False), (True, False), (True, True)]:
    for normal in range(3):
        for positive in [False, True]:
            full = run(normal, positive, False, refined, graded)
            half = run(normal, positive, True, refined, graded)
            errors = {}
            for field in ['e', 'h']:
                # The compact H(curl) reconstruction supplies tangential E
                # at the weak PMC face. Normal E and tangential H there are
                # unilateral cell samples, not symmetric trace averages.
                # Compare the well-defined tangential E trace and full E/H
                # vectors at all three matched interior points.
                indexes = range(4) if field == 'e' else range(1, 4)
                peak = max(float(np.max(abs(full[field+str(i)][:, 1:]))) for i in indexes)
                assert peak > 1e-12, (normal, refined, field, 'unexcited comparison')
                error = 0.
                for i in indexes:
                    key = field+str(i)
                    np.testing.assert_allclose(full[key][:,0], half[key][:,0], rtol=1e-12, atol=0)
                    assert np.isfinite(half[key]).all(), key
                    columns = [n+1 for n in range(3) if n != normal] if field == 'e' and i == 0 else [1,2,3]
                    error = max(error, float(np.max(abs(full[key][:,columns]-half[key][:,columns])))/peak)
                errors[field] = error
                assert error < 2e-4, (normal, positive, refined, field, error)
            label = '{}{}-{}'.format('xyz'[normal], '+' if positive else '-', ('local-graded' if graded else 'local') if refined else 'identity')
            metrics[label] = errors
            print('PMC mirrored-domain comparison', label, errors, flush=True)

# Adjacent Mur corners and PML must still fail closed.
for condition in ['MUR', 'PML_8']:
    fdtd = openEMS(NrTS=1, CartesianSubgridBox=[-4,-4,-4,4,4,4])
    fdtd.SetBoundaryCond([condition]*6)
    fdtd.SetGaussExcite(13e9, 9e9)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    for axis in 'xyz':
        csx.GetGrid().SetLines(axis, np.arange(-8,9,2))
    csx.GetGrid().SetDeltaUnit(1e-3)
    old = os.getcwd()
    try:
        try:
            fdtd.Run(os.path.join(Sim_Path, 'rejected-'+condition), engine='basic', setup_only=True)
        except (ValueError, RuntimeError) as error:
            assert 'CartesianSubgrid' in str(error), str(error)
        else:
            raise AssertionError('unsupported boundary accepted: '+condition)
    finally:
        os.chdir(old)
with open(os.path.join(Sim_Path, 'metrics.json'), 'w') as f:
    json.dump(metrics, f, indent=2)
print('All PMC native symmetry controls passed:', metrics, flush=True)

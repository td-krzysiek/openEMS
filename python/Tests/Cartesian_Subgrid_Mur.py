# SPDX-License-Identifier: GPL-3.0-or-later
"""Parallel Mur faces: Yee equivalence and finite-window TE reflection.

The long PEC reference has no return in the measurement window. At oblique
incidence a fixed transverse wave number is used: the quoted angle is at
10 GHz, not constant over frequency. This is a first-order ABC, not PML.
"""
import json
import os
import tempfile

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_Mur')
os.makedirs(Sim_Path, exist_ok=True)
C0 = 299792458.
metrics = {}


def run(angle, compact, long=False, refined=False, graded=False):
    suffix = ('-local' if refined else '') + ('-graded' if graded else '')
    label = '{}-{}-{}{}'.format(angle, 'compact' if compact else 'yee',
                               'long' if long else 'mur', suffix)
    fdtd = openEMS(NrTS=12000, EndCriteria=0, TimeStep=5e-13, OverSampling=1)
    # Keep illumination away from the oblique mode cutoff; a wider pulse
    # has a long near-cutoff tail that failed the finite-window check.
    fdtd.SetGaussExcite(10e9, 1e9)
    bc = ['PEC'] * 6
    bc[0] = 'MUR'
    if not long:
        bc[1] = 'MUR'
    if angle == 0:
        bc[4] = bc[5] = 'PMC'
    fdtd.SetBoundaryCond(bc)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    width = 4. if angle == 0 else C0 / (2 * 10e9 * np.sin(np.deg2rad(angle))) * 1e3
    x = np.arange(0, 1001 if long else 121, 1.)
    if graded:
        x[1::2] -= .2
    grid.SetLines('x', x)
    grid.SetLines('y', np.linspace(0, width, 5 if angle == 0 else 17))
    grid.SetLines('z', np.arange(5.))
    if compact:
        fdtd.AddCartesianSubgrid(box=[50, width/4, 1, 80, 3*width/4, 3],
                                ratio=[2 if refined else 1, 1, 1])
    vector = [0, 1, 0] if angle == 0 else [0, 0, 1]
    source = csx.AddExcitation('source', exc_type=0, exc_val=vector)
    if angle:
        source.SetWeightFunction(['1', '1', 'sin(pi*y/{:.17g})'.format(width)])
    source.AddBox([30, 0, 0], [30, width, 4])
    csx.AddProbe('field', p_type=2).AddPoint([40, width/2, 2])
    path = os.path.join(Sim_Path, label)
    os.makedirs(path, exist_ok=True)
    fdtd.Write2XML(os.path.join(path, 'input.xml'))
    old = os.getcwd()
    try:
        fdtd.Run(path, engine='basic', numThreads=1, exact_endcriteria=True)
    finally:
        os.chdir(old)
    data = np.loadtxt(os.path.join(path, 'field'), comments='%')
    assert np.isfinite(data).all(), label
    return data[:, 0], data[:, 2 if angle == 0 else 3]


for angle, graded in [(0, False), (30, False), (45, False), (30, True)]:
    t, reference = run(angle, False, long=True, graded=graded)
    _, yee = run(angle, False, graded=graded)
    _, compact = run(angle, True, graded=graded)
    peak = float(np.max(abs(reference)))
    assert peak > 1e-12, 'unexcited reference'
    error = float(np.max(abs(compact-yee))/peak)
    assert error < 2e-4, (angle, 'Yee equivalence', error)
    # The matched left Mur face preserves source-side returns. The remote
    # PEC end cannot return before the 6 ns observation cutoff. Retain the
    # complete incident pulse rather than truncating its late oscillations.
    incident = reference
    reflected = compact-reference
    frequencies = np.array([9.5e9, 10e9, 10.5e9])
    kernel = np.exp(-2j*np.pi*frequencies[:, None]*t)
    inc = kernel @ incident
    ref = kernel @ reflected
    assert np.min(abs(inc)) > .1*np.max(abs(inc)), 'poor spectral illumination'
    reflection = abs(ref/inc)
    cosine = np.sqrt(1-(10e9*np.sin(np.deg2rad(angle))/frequencies)**2)
    continuum = (1-cosine)/(1+cosine)
    shorter = t < 4e-9
    short_reflection = abs((kernel[:, shorter] @ reflected[shorter]) /
                           (kernel[:, shorter] @ incident[shorter]))
    window_change = float(np.max(abs(reflection-short_reflection)))
    assert window_change < .01, (angle, 'unsettled spectral window', window_change)
    label = str(angle) + ('-graded' if graded else '')
    metrics[label] = dict(yee_peak_error=error,
                          window_change_4_to_6_ns=window_change,
                          frequencies_hz=frequencies.tolist(),
                          continuum_first_order_magnitude=continuum.tolist(),
                          reflection_magnitude=reflection.tolist())
    # Mur is deliberately approximate at oblique incidence; these gates
    # catch broken absorption while preserving measured angle dependence.
    assert np.max(reflection) < (.08 if angle == 0 else .4), (angle, reflection)
    print('Mur reflection', angle, metrics[label], flush=True)

# Local refinement: the matched long compact reference subtracts the same
# coarse/fine scattering. This separates the end-boundary reflection.
t, reference = run(0, True, long=True, refined=True)
_, local = run(0, True, refined=True)
peak = np.max(abs(reference))
returned = float(np.max(abs((local-reference)[t >= .5e-9]))/peak)
assert returned < .08, returned
metrics['local'] = dict(returned_peak=returned)
with open(os.path.join(Sim_Path, 'metrics.json'), 'w') as stream:
    json.dump(metrics, stream, indent=2)
print('All native Mur controls passed', metrics, flush=True)


def reject(label, alter, reason):
    fdtd = openEMS(NrTS=1, CartesianSubgridBox=[3, 3, 3, 5, 5, 5])
    fdtd.SetGaussExcite(10e9, 1e9)
    fdtd.SetBoundaryCond(['MUR', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC'])
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    csx.GetGrid().SetDeltaUnit(1e-3)
    for axis in 'xyz':
        csx.GetGrid().SetLines(axis, np.arange(9.))
    alter(fdtd, csx)
    path = os.path.join(Sim_Path, 'reject-'+label)
    os.makedirs(path, exist_ok=True)
    fdtd.Write2XML(os.path.join(path, 'input.xml'))
    old = os.getcwd()
    try:
        try:
            fdtd.Run(path, engine='basic', setup_only=True)
        except (ValueError, RuntimeError) as error:
            assert reason in str(error), (label, str(error), reason)
        else:
            raise AssertionError('unsupported Mur configuration accepted: '+label)
    finally:
        os.chdir(old)


reject('adjacent', lambda f, c: f.SetBoundaryCond(['MUR', 'PEC', 'MUR', 'PEC', 'PEC', 'PEC']),
       'parallel')
reject('refined-collar', lambda f, c: f.SetCartesianSubgrids([
    dict(box=[1, 3, 3, 3, 5, 5], ratio=[2, 2, 2])]), 'two unrefined')
reject('heterogeneous-collar', lambda f, c: c.AddMaterial('slab', epsilon=2).AddBox(
    [0, 0, 0], [1, 8, 8]), 'homogeneous')
reject('metal-collar', lambda f, c: c.AddMetal('metal').AddBox([0, 3, 3], [1, 5, 5]),
       'PEC touches')
reject('source-collar', lambda f, c: c.AddExcitation('source', exc_type=0,
       exc_val=[0, 0, 1]).AddBox([1, 3, 3], [1, 3, 4]), 'outside the Mur collar')
reject('resistor-collar', lambda f, c: c.AddLumpedElement('R', ny=2, R=50).AddBox(
    [1, 3, 3], [1, 3, 4]), 'outside the Mur collar')
print('Mur scope guards passed', flush=True)

# SPDX-License-Identifier: GPL-3.0-or-later
"""Classical convolution PML: finite-window reflection and scope rejection.

The long reference uses the identical compact scheme/source-side PML, with
the far boundary outside the measurement light cone. This isolates boundary
returns from coarse/fine scattering. Angles refer to 10 GHz only.
"""
import json
import os
import tempfile

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_PML')
os.makedirs(Sim_Path, exist_ok=True)
C0 = 299792458.
metrics = {}


def run(angle=0, polarization='TE', long=False, refined=False, graded=False, pec=False, yee=False):
    label = '{}-{}-{}-{}-{}-{}-{}'.format(angle, polarization, 'long' if long else 'short',
        'local' if refined else 'identity', 'graded' if graded else 'uniform',
        'pec' if pec else 'pml', 'yee' if yee else 'compact')
    fdtd = openEMS(NrTS=12000, EndCriteria=1e-30, TimeStep=5e-13, OverSampling=1)
    fdtd.SetGaussExcite(10e9, 1e9)
    bc = ['PEC'] * 6
    bc[0] = 'PML_8'
    if not long and not pec:
        bc[1] = 'PML_8'
    if angle == 0:
        bc[4] = bc[5] = 'PMC'
    elif polarization == 'TM':
        bc[2] = bc[3] = 'PMC'
        bc[4] = bc[5] = 'PMC'
    fdtd.SetBoundaryCond(bc)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    width = 4. if angle == 0 else C0/(2*10e9*np.sin(np.deg2rad(angle)))*1e3
    x = np.arange(0, 1001 if long else 121, 1.)
    if graded:
        x[1::2] -= .2
    mesh.SetLines('x', x)
    mesh.SetLines('y', np.linspace(0, width, 5 if angle == 0 else 17))
    mesh.SetLines('z', np.arange(5.))
    if not yee:
        fdtd.AddCartesianSubgrid(box=[50, width/4, 1, 80, 3*width/4, 3],
            ratio=[2 if refined else 1, 1, 1])
    # TE uses Ez with PEC side walls; TM uses Ey with PMC side walls.
    # Both are tangential to the x=30 source plane and resolve actual edges.
    vector = [0, 0, 1] if angle and polarization == 'TE' else [0, 1, 0]
    source = csx.AddExcitation('source', exc_type=0, exc_val=vector)
    if angle:
        weights = ['1']*3
        weights[2 if polarization == 'TE' else 1] = 'sin(pi*y/{:.17g})'.format(width)
        source.SetWeightFunction(weights)
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
    # Probes sample at the excitation Nyquist interval, not every FDTD step.
    sample_interval = data[1, 0]-data[0, 0]
    assert len(data) > 100 and 6e-9-sample_interval-1e-16 <= data[-1, 0] <= 6e-9+1e-16, (
        label, 'record does not reach the fixed physical duration', data[-1, 0])
    return data[:, 0], data[:, 3 if angle and polarization == 'TE' else 2]


for angle, polarization, graded, refined in [
        (0, 'TE', False, False), (30, 'TE', False, False),
        (45, 'TE', False, False), (30, 'TM', False, False),
        (30, 'TE', True, False), (0, 'TE', False, True)]:
    t, incident = run(angle, polarization, long=True, graded=graded, refined=refined)
    _, short = run(angle, polarization, graded=graded, refined=refined)
    reflected = short-incident
    frequencies = np.array([9.5e9, 10e9, 10.5e9])
    kernel = np.exp(-2j*np.pi*frequencies[:, None]*t)
    inc = kernel @ incident
    assert np.max(abs(incident)) > 1e-12, 'unexcited reference'
    assert np.min(abs(inc)) > .1*np.max(abs(inc)), 'poor spectral illumination'
    reflection = abs((kernel @ reflected)/inc)
    earlier = t < 4e-9
    short_window = abs((kernel[:, earlier] @ reflected[earlier]) /
                       (kernel[:, earlier] @ incident[earlier]))
    window_change = float(np.max(abs(reflection-short_window)))
    label = '{}-{}-graded{}-local{}'.format(angle, polarization, graded, refined)
    metrics[label] = dict(frequencies_hz=frequencies.tolist(), reflection_magnitude=reflection.tolist(),
        window_change_4_to_6_ns=window_change, returned_peak=float(np.max(abs(reflected))/np.max(abs(incident))))
    print('PML reflection', label, metrics[label], flush=True)
    # Conservative qualification gates: <=2% amplitude boundary reflection
    # and <=1 percentage point finite-window change across the sampled band.
    assert window_change < .01, (label, 'unsettled window', window_change)
    assert np.max(reflection) < .02, (label, 'boundary reflection', reflection)
    if angle == 0 and not refined:
        _, negative = run(pec=True)
        negative_reflection = abs((kernel @ (negative-incident))/inc)
        assert np.min(negative_reflection) > .8, ('negative PEC control', negative_reflection)
        metrics['negative-pec'] = negative_reflection.tolist()
        # Ordinary Yee agrees in the early incident window, before waves
        # can reach either source-side absorber and return to the probe.
        _, ordinary = run(yee=True)
        before_return = t < .7e-9
        error = float(np.max(abs(ordinary[before_return]-short[before_return]))/np.max(abs(incident)))
        assert error < 2e-4, ('incident Yee comparison', error)
        metrics['incident-yee-peak-error'] = error


def reject(label, alter, reason):
    fdtd = openEMS(NrTS=1, CartesianSubgridBox=[8, 8, 8, 12, 12, 12])
    fdtd.SetGaussExcite(10e9, 1e9)
    fdtd.SetBoundaryCond(['PML_4']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    csx.GetGrid().SetDeltaUnit(1e-3)
    for axis in 'xyz':
        csx.GetGrid().SetLines(axis, np.arange(21.))
    alter(fdtd, csx)
    path = os.path.join(Sim_Path, 'reject-'+label)
    os.makedirs(path, exist_ok=True)
    fdtd.Write2XML(os.path.join(path, 'input.xml'))
    old = os.getcwd()
    try:
        try:
            fdtd.Run(path, engine='basic', setup_only=True)
        except (RuntimeError, ValueError) as error:
            assert reason in str(error), (label, str(error), reason)
        else:
            raise AssertionError('unsupported PML configuration accepted: '+label)
    finally:
        os.chdir(old)


reject('thin', lambda f,c: f.SetBoundaryCond(['PML_0']*6), 'four base cells')
reject('overlap', lambda f,c: f.SetBoundaryCond(['PML_10']*6), 'intervening')
reject('mur-mix', lambda f,c: f.SetBoundaryCond(['MUR','PML_4','PEC','PEC','PEC','PEC']), 'mixed Mur')
reject('refinement-collar', lambda f,c: f.SetCartesianSubgrids([
    dict(box=[4,8,8,6,12,12], ratio=[2,2,2])]), 'one-cell guard')
reject('material-collar', lambda f,c: c.AddMaterial('slab', epsilon=2).AddBox([0,0,0],[1,20,20]), 'homogeneous')
reject('pec-collar', lambda f,c: c.AddMetal('metal').AddBox([1,8,8],[2,12,12]), 'internal PEC')
reject('source-collar', lambda f,c: c.AddExcitation('drive', exc_type=0, exc_val=[0,1,0]).AddBox(
    [2,8,8],[2,12,12]), 'sources and resistors')
reject('resistor-collar', lambda f,c: c.AddLumpedElement('load', ny=1, R=50).AddBox(
    [2,8,8],[2,9,8]), 'sources and resistors')

with open(os.path.join(Sim_Path, 'metrics.json'), 'w') as stream:
    json.dump(metrics, stream, indent=2)
print('All native PML controls passed', metrics, flush=True)

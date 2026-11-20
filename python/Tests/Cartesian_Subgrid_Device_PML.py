# SPDX-License-Identifier: GPL-3.0-or-later
"""Matched-input PML CPU/device E/H, FD, graded/local and TE/TM controls.

The same small fixtures as Cartesian_Subgrid_PML.py, with deterministic XML
and extra E/H probes/dumps. The original CPU reflection/window controls remain
separate and unchanged. This is not yet a device antenna/NF2FF qualification.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from _subgrid_device_support import run_native_device

SIM = Path(tempfile.gettempdir()) / 'Cartesian_Subgrid_Device_PML'
SIM.mkdir(parents=True, exist_ok=True)
PHYSICAL = os.environ.get('OPENEMS_SUBGRID_DEVICE', '')
assert PHYSICAL in ('', 'metal', 'cuda')
C0 = 299792458.
GATE = 2e-4


def fixture(angle, polarization, graded, refined):
    fdtd = openEMS(NrTS=12000, EndCriteria=0, TimeStep=5e-13, OverSampling=1)
    fdtd.SetGaussExcite(10e9, 1e9)
    bc = ['PML_8', 'PML_8', 'PEC', 'PEC', 'PEC', 'PEC']
    if angle == 0:
        bc[4] = bc[5] = 'PMC'
    elif polarization == 'TM':
        bc[2:] = ['PMC'] * 4
    fdtd.SetBoundaryCond(bc)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    width = 4. if angle == 0 else C0 / (2 * 10e9 * np.sin(np.deg2rad(angle))) * 1e3
    x = np.arange(121.)
    if graded:
        x[1::2] -= .2
    mesh.SetLines('x', x)
    mesh.SetLines('y', np.linspace(0, width, 5 if angle == 0 else 17))
    mesh.SetLines('z', np.arange(5.))
    fdtd.AddCartesianSubgrid(box=[50, width / 4, 1, 80, 3 * width / 4, 3],
                            ratio=[2 if refined else 1, 1, 1])
    vector = [0, 0, 1] if angle and polarization == 'TE' else [0, 1, 0]
    source = csx.AddExcitation('source', exc_type=0, exc_val=vector)
    if angle:
        weights = ['1'] * 3
        weights[2 if polarization == 'TE' else 1] = f'sin(pi*y/{width:.17g})'
        source.SetWeightFunction(weights)
    source.AddBox([30, 0, 0], [30, width, 4])
    for name, quantity in [('electric', 2), ('magnetic', 3)]:
        csx.AddProbe(name, p_type=quantity).AddPoint([40, width / 2, 2])
    for name, quantity in [('ef', 10), ('hf', 11)]:
        csx.AddDump(name, dump_type=quantity, file_type=1, frequency=[9.5e9, 10e9, 10.5e9]).AddBox(
            [40, 0, 0], [40, width, 4])
    for prop in csx.GetAllProperties():
        prop.SetColor([96, 128, 160])
    return fdtd


def run(label, engine, options):
    fdtd = fixture(*options)
    path = SIM / (label + '-' + engine)
    path.mkdir(exist_ok=False)
    fdtd.Write2XML(str(path / 'input.xml'))
    previous = os.getcwd()
    try:
        result = run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1, exact_endcriteria=True)
        assert result in (None, 0), (label, engine, result)
    finally:
        os.chdir(previous)
    return path


def compare(reference, candidate):
    assert hashlib.sha256((reference / 'input.xml').read_bytes()).digest() == hashlib.sha256(
        (candidate / 'input.xml').read_bytes()).digest(), 'different native PML inputs'
    errors = {}
    for name in ['electric', 'magnetic']:
        a, b = [np.loadtxt(path / name, comments='%') for path in [reference, candidate]]
        assert a.shape == b.shape and len(a) > 100 and np.isfinite(b).all()
        np.testing.assert_allclose(a[:, 0], b[:, 0], rtol=1e-12, atol=0)
        scale = float(np.max(abs(a[:, 1:])))
        assert scale > 1e-12
        errors[name] = float(np.max(abs(a[:, 1:] - b[:, 1:]))) / scale
    for name in ['ef.h5', 'hf.h5']:
        with h5py.File(reference / name) as a, h5py.File(candidate / name) as b:
            keys = []
            a['FieldData'].visititems(lambda key, obj: keys.append(key) if isinstance(obj, h5py.Dataset) else None)
            assert keys
            scale = max(float(np.max(abs(a['FieldData/' + key][...]))) for key in keys)
            delta = 0
            for key in keys:
                left, right = a['FieldData/' + key][...], b['FieldData/' + key][...]
                assert left.shape == right.shape and np.isfinite(right).all()
                delta = max(delta, float(np.max(abs(left - right))))
            assert scale > 1e-30
            errors[name] = delta / scale
    assert max(errors.values()) < GATE, (candidate.name, errors)
    print(candidate.name, errors, flush=True)
    return errors


metrics = {}
for options in [(0, 'TE', False, False), (30, 'TE', False, False), (45, 'TE', False, False),
                (30, 'TM', False, False), (30, 'TE', True, False), (0, 'TE', False, True)]:
    label = '{}-{}-graded{}-local{}'.format(*options)
    cpu = run(label, 'basic', options)
    reference = run(label, 'gpu-reference', options)
    metrics[label + '-reference'] = compare(cpu, reference)
    if PHYSICAL:
        device = run(label, 'gpu', options)
        metrics[label + '-' + PHYSICAL] = compare(cpu, device)
        metrics[label + '-' + PHYSICAL + '-v-reference'] = compare(reference, device)
    (SIM / 'metrics.json').write_text(json.dumps({'physical': PHYSICAL, 'gate': GATE, 'errors': metrics}, indent=2) + '\n')
print('PASS: native paired device PML E/H, frequency dumps, TE/TM, local and graded controls', flush=True)

# SPDX-License-Identifier: GPL-3.0-or-later
"""Native PEC polygons, finite Cu18 contacts and strict axial alignment.

Uses CSXCAD's existing planar/extruded polygon membership, without changing the
continuous contours into rectangle unions. The grid still samples the contour;
this is neither a conformal-boundary method nor fabrication/RF convergence.
A native point probe interpolates staggered base samples, so PEC-exclusion
witnesses below keep its entire interpolation footprint inside copper.
"""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from _subgrid_device_support import run_native_device

ROOT = Path(tempfile.gettempdir()) / 'Cartesian_Subgrid_Polygons'
ROOT.mkdir(exist_ok=False)
PHYSICAL = os.environ.get('OPENEMS_SUBGRID_DEVICE', '')
assert PHYSICAL in ('', 'metal', 'cuda')
ENGINES = ['basic', 'gpu-reference'] + (['gpu'] if PHYSICAL else [])
GATE = 2e-4


def contour(shape, count=64):
    if shape in ('rectangle', 'box', 'negative', 'flat', 'flat-box'):
        return np.array([[-2, 2, 2, -2], [-2, -2, 2, 2]], dtype=float)
    if shape == 'triangle':
        return np.array([[-2, 2, 0], [-2, -2, 2]], dtype=float)
    if shape in ('disk', 'cylinder'):
        angle = np.linspace(0, 2*np.pi, count, endpoint=False)
        return 2.47*np.array([np.cos(angle), np.sin(angle)])
    if shape == 'bend':
        outer = np.linspace(np.pi, 0, count+1)
        inner = np.linspace(0, np.pi, count+1)
        return np.column_stack(([-3, -3], 3*np.array([np.cos(outer), np.sin(outer)]),
                                [3, -3], [1, -3], np.array([np.cos(inner), np.sin(inner)]),
                                [-1, -3]))
    raise ValueError(shape)


def fixture(shape, normal=2, offset=(0, 0, 0), scale=1, reverse=False, thin=False, count=64):
    shift = (normal-2) % 3
    offset = np.asarray(offset, dtype=float)

    def point(value):
        return ((np.roll(np.asarray(value, dtype=float), shift)+offset)*scale).tolist()

    steps = 20000 if thin else 2048
    fdtd = openEMS(NrTS=steps, EndCriteria=0, TimeStepFactor=.8)
    fdtd.SetCartesianSubgrid(point([-4, -4, -2 if thin else -4]) + point([4, 4, 2 if thin else 4]))
    # Corners are ordered by coordinate; cyclic rotation preserves low/high.
    fdtd.SetCartesianSubgridRatio(np.roll([2, 2, 1] if thin else [2, 2, 2], shift).tolist())
    fdtd.SetGaussExcite(13e9, 9e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    csx.GetGrid().SetDeltaUnit(1e-3/scale)
    for axis in range(3):
        old_axis = (axis-shift) % 3
        lines = [-8, -4, -2, -.25, 0, .018, 2, 4, 8] if thin and old_axis == 2 else np.arange(-8, 9, 2)
        csx.GetGrid().SetLines('xyz'[axis], ((np.asarray(lines)+offset[axis])*scale).tolist())
    direction = np.roll([0, 0, 1], shift).tolist()
    csx.AddExcitation('source', exc_type=0, exc_val=direction).AddBox(point([-6, 0, -2]), point([-6, 0, 2]))
    csx.AddProbe('received', p_type=2).AddPoint(point([6, 4, 2]))
    csx.AddProbe('magnetic', p_type=3).AddPoint(point([6, 4, 2]))
    inside = [-2, -2, .009 if thin else 0] if shape == 'bend' else [0, 0, .009 if thin else 0]
    csx.AddProbe('local_e', p_type=2).AddPoint(point(inside))
    csx.AddProbe('local_h', p_type=3).AddPoint(point(inside))
    csx.AddMaterial('dielectric', epsilon=2.25).AddBox(
        point([-8, -8, -.25 if thin else -8]), point([8, 8, 0]), priority=1)
    metal = csx.AddMetal('conductor')
    lo, hi = (0, .018) if thin else (-3, 3)
    if shape in ('flat', 'flat-box'):
        lo = hi = 0
    vertices = contour(shape, count)
    if reverse:
        vertices = vertices[:, ::-1]
    plane_axes = [(normal+1) % 3, (normal+2) % 3]
    vertices = (vertices+offset[plane_axes, None])*scale
    if shape in ('box', 'flat-box'):
        metal.AddBox(point([-2, -2, lo]), point([2, 2, hi]), priority=10)
    elif shape == 'cylinder':
        metal.AddCylinder(point([0, 0, lo]), point([0, 0, hi]), 2.47*scale, priority=10)
    elif shape == 'flat':
        metal.AddPolygon(points=vertices, norm_dir=normal, elevation=(lo+offset[normal])*scale, priority=10)
    else:
        elevation, length = ((hi, lo-hi) if shape == 'negative' else (lo, hi-lo))
        metal.AddLinPoly(points=vertices, norm_dir=normal,
                         elevation=(elevation+offset[normal])*scale, length=length*scale, priority=10)
    for prop in csx.GetAllProperties():
        prop.SetColor([96, 128, 160])
    return fdtd, csx, shift


def data(path):
    result = {}
    for name in ('received', 'magnetic', 'local_e', 'local_h'):
        values = np.loadtxt(path/name, comments='%')
        assert values.ndim == 2 and values.shape[0] > 10 and np.isfinite(values).all(), (path, name)
        result[name] = values
    return result


def compare(a, b, gate, label, shift=0):
    errors = {}
    for key in a:
        reference, candidate = a[key], b[key].copy()
        np.testing.assert_allclose(reference[:, 0], candidate[:, 0], rtol=1e-11, atol=0)
        if shift:
            candidate[:, 1:] = np.roll(candidate[:, 1:], -shift, axis=1)
        peak = max(np.max(abs(reference[:, 1:])), np.max(abs(candidate[:, 1:])))
        delta = np.max(abs(reference[:, 1:]-candidate[:, 1:]))
        error = delta / max(peak, 1e-300)
        assert error < gate or delta < 1e-13, (label, key, error, delta)
        errors[key] = float(error)
    return errors


metrics, results = {}, {}
cases = [('box', {}), ('rectangle', {}), ('reverse', {'shape': 'rectangle', 'reverse': True}),
         ('negative', {}), ('triangle', {}), ('cylinder', {}), ('disk', {}),
         ('bend', {}), ('translated', {'shape': 'bend', 'offset': [31, -23, 17]}),
         ('micrometre', {'shape': 'bend', 'scale': 1000}),
         ('normal-x', {'shape': 'rectangle', 'normal': 0}),
         ('normal-y', {'shape': 'rectangle', 'normal': 1}),
         ('flat-box', {}), ('flat', {}), ('cu18-bend', {'shape': 'bend', 'thin': True})]
for label, options in cases:
    for engine in ENGINES:
        config = dict(shape=label)
        config.update(options)
        path = ROOT/(label+'-'+engine)
        path.mkdir(exist_ok=False)
        fdtd, csx, shift = fixture(**config)
        fdtd.Write2XML(str(path/'input.xml'))
        previous = os.getcwd()
        try:
            run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1, exact_endcriteria=True)
        finally:
            os.chdir(previous)
        current = data(path)
        assert np.max(abs(current['received'][:, 1:])) > 1e-7, (label, engine, 'no exterior response')
        results[label, engine] = current
        row = {'inputSha256': hashlib.sha256((path/'input.xml').read_bytes()).hexdigest()}
        if engine != 'basic':
            reference = ROOT/(label+'-basic')
            assert (path/'input.xml').read_bytes() == (reference/'input.xml').read_bytes(), label
            row['errors'] = compare(results[label, 'basic'], current, GATE, label+'-'+engine)
        metrics[label+'-'+engine] = row
        print(label, engine, row, flush=True)

for label, parent in [('rectangle', 'box'), ('reverse', 'box'), ('negative', 'box'),
                      ('disk', 'cylinder'), ('flat', 'flat-box'), ('translated', 'bend'), ('micrometre', 'bend')]:
    metrics[label+'-equivalence'] = compare(results[parent, 'basic'], results[label, 'basic'], 1e-8, label)
for label, shift in [('normal-x', 1), ('normal-y', 2)]:
    metrics[label+'-rotation'] = compare(results['rectangle', 'basic'], results[label, 'basic'], 1e-8, label, shift)
(ROOT/'metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
print('PASS: native finite PEC polygon fields and equivalent shapes/axes/units', flush=True)


# Strict material/axial alignment: coordinate roundoff is not geometry snapping.
for name, options in [('box', {'shape': 'box'}), ('bend', {'shape': 'bend'}),
                      ('negative', {'shape': 'negative'}), ('cylinder', {'shape': 'cylinder'}),
                      ('translated', {'shape': 'bend', 'offset': [31, -23, 17]}),
                      ('units', {'shape': 'bend', 'scale': 1000}),
                      ('normal-x', {'shape': 'bend', 'normal': 0}),
                      ('normal-y', {'shape': 'bend', 'normal': 1})]:
    fdtd, csx, shift = fixture(thin=True, **options)
    path = ROOT / ('aligned-' + name)
    fdtd.Write2XML(str(path) + '.xml')
    fdtd.Run(str(path), engine='basic', setup_only=True, numThreads=1)
    print('accepted aligned', name, flush=True)

for sign in [-1, 1]:
    for scale in [1, 1000]:
        fdtd, csx, shift = fixture('box', thin=True, scale=scale)
        conductor = next(prop for prop in csx.GetAllProperties() if prop.GetName() == 'conductor')
        conductor.GetPrimitive(0).SetStop([2*scale, 2*scale, (.018 + sign*1e-9)*scale])
        path = ROOT / f'cut-{sign}-{scale}'
        fdtd.Write2XML(str(path) + '.xml')
        try:
            fdtd.Run(str(path), engine='basic', setup_only=True, numThreads=1)
        except (RuntimeError, ValueError) as error:
            assert 'boundary cuts a leaf' in str(error), error
            print('rejected actual boundary', sign, scale, error, flush=True)
        else:
            raise AssertionError('A real 1e-9 mm displacement was silently snapped')
print('PASS: thin axial alignment, axes/units/origin and strict negative controls')


# Exclusion probes: resolved interpolation support, including the Cu18 interior.
exclusion_metrics = {}
for thin in [False, True]:
    for engine in ENGINES:
        fdtd, csx, shift = fixture('bend', thin=thin)
        for axis in 'xy':
            csx.GetGrid().SetLines(axis, [-8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8])
        if thin:
            csx.GetGrid().SetLines('z', [-8, -4, -2, -.25, 0, .009, .018, 2, 4, 8])
        path = ROOT / f'exclusion-{thin}-{engine}'
        path.mkdir()
        fdtd.Write2XML(str(path / 'input.xml'))
        run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1, exact_endcriteria=True)
        row = {}
        for name in ['local_e', 'local_h']:
            values = np.loadtxt(path/name, comments='%')
            peak = float(np.max(abs(values[:, 1:])))
            assert peak < 1e-12, (thin, engine, name, peak)
            row[name] = peak
        exterior = np.loadtxt(path/'received', comments='%')
        assert np.max(abs(exterior[:, 1:])) > 1e-7
        exclusion_metrics[f'{thin}-{engine}'] = row
        print(thin, engine, row, flush=True)
(ROOT/'exclusion-metrics.json').write_text(json.dumps(exclusion_metrics, indent=2)+'\n')
print('PASS: resolved native-probe footprint excludes fields from solid/18um U bends')


# A finite-height strip and ground retain the existing one-edge port contract.
# This fixture checks geometry/contact mapping, not production 24GHz calibration.
def port_fixture(shape):
    fdtd = openEMS(NrTS=24000, EndCriteria=0, TimeStepFactor=1,
                  TimeStep=5e-14, OverSampling=20)
    fdtd.SetCartesianSubgrid([1.2, -.8, .2, 2.8, .8, 1.0])
    fdtd.SetCartesianSubgridRatio([2, 1, 1])
    fdtd.SetGaussExcite(24e9, 12e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis, limits in [('x', (0, 4)), ('y', (-1.2, 1.2))]:
        grid.SetLines(axis, np.round(np.linspace(*limits, int(round((limits[1]-limits[0])/.1))+1), 12))
    grid.SetLines('z', [0, .182, .2, .4, .418, .6, .8, 1, 1.2])
    csx.AddMaterial('substrate', epsilon=2.25).AddBox([0, -1.2, .2], [4, 1.2, .4])
    metal = csx.AddMetal('electrodes_cu18')
    metal.AddBox([0, -1.2, .182], [4, 1.2, .2], priority=10)
    if shape == 'box':
        metal.AddBox([.6, -.2, .4], [3.4, .2, .418], priority=10)
    else:
        metal.AddLinPoly(points=[[.6, 3.4, 3.4, .6], [-.2, -.2, .2, .2]],
                         norm_dir=2, elevation=.4, length=.018, priority=10)
    ports = [fdtd.AddLumpedPort(number, 50, [x, 0, .2], [x, 0, .4], 'z',
                               1 if number == 1 else 0)
             for number, x in enumerate([.6, 3.4], 1)]
    for prop in csx.GetAllProperties():
        prop.SetColor([96, 128, 160])
    return fdtd, csx, ports


port_metrics = {}
port_reference = None
port_reference_values = None
for shape in ['box', 'polygon']:
    for engine in ENGINES:
        fdtd, csx, ports = port_fixture(shape)
        path = ROOT / f'ports-{shape}-{engine}'
        path.mkdir()
        fdtd.Write2XML(str(path / 'input.xml'))
        if engine != 'basic':
            assert (path/'input.xml').read_bytes() == (ROOT/f'ports-{shape}-basic'/'input.xml').read_bytes()
        run_native_device(fdtd, path, engine, PHYSICAL, numThreads=1, exact_endcriteria=True)
        waves = {name: np.loadtxt(path/name, comments='%')
                 for name in ['port_ut_1', 'port_it_1', 'port_ut_2', 'port_it_2']}
        assert all(np.isfinite(values).all() for values in waves.values())
        assert np.max(abs(waves['port_ut_1'][:, 1])) > 1e-7
        if port_reference is None:
            port_reference = waves
        errors = {}
        for name, values in waves.items():
            np.testing.assert_allclose(values[:, 0], port_reference[name][:, 0], rtol=1e-11, atol=0)
            error = np.max(abs(values[:, 1]-port_reference[name][:, 1])) / np.max(abs(port_reference[name][:, 1]))
            assert error < GATE, (shape, engine, name, error)
            errors[name] = float(error)
        for port in ports:
            port.CalcPort(str(path), np.linspace(20e9, 28e9, 41))
        values = dict(z=ports[0].uf_tot/ports[0].if_tot, s11=ports[0].uf_ref/ports[0].uf_inc,
                      s21=ports[1].uf_ref/ports[0].uf_inc)
        assert all(np.isfinite(value).all() for value in values.values())
        if port_reference_values is None:
            port_reference_values = values
        for name, value in values.items():
            error = float(np.max(abs(value-port_reference_values[name])))
            assert error < (.01 if name == 'z' else GATE), (shape, engine, name, error)
            errors[name] = error
        phase = float(np.max(abs(np.angle(values['s21']/port_reference_values['s21'], deg=True))))
        assert phase < .01, (shape, engine, phase)
        errors['phase_deg'] = phase
        balance = np.max(abs(abs(values['s11'])**2 + abs(values['s21'])**2 - 1))
        assert balance < 1e-3, (shape, engine, balance)
        port_metrics[f'{shape}-{engine}'] = {'wave_errors': errors, 'power_error': float(balance)}
        print(shape, engine, port_metrics[f'{shape}-{engine}'], flush=True)
(ROOT/'port-metrics.json').write_text(json.dumps(port_metrics, indent=2)+'\n')
print('PASS: native PEC polygon/Cu18 fields, contacts, power, axes/units and alignment')


def reject(label, change, expected):
    fdtd, csx, shift = fixture('box')
    change(csx)
    path = ROOT / ('guard-' + label)
    fdtd.Write2XML(str(path) + '.xml')
    try:
        fdtd.Run(str(path), engine='basic', setup_only=True, numThreads=1)
    except (RuntimeError, ValueError) as error:
        assert 'CartesianSubgrid' in str(error) and expected in str(error), (label, error)
        print('rejected', label, error, flush=True)
        return
    raise AssertionError('Invalid polygon setup accepted: ' + label)


square = [[-2, 2, 2, -2], [-2, -2, 2, 2]]
reject('too-few', lambda c: c.AddMetal('bad').AddPolygon(
    points=[[-2, 2], [-2, 2]], norm_dir=2, elevation=0), 'at least 3 vertices')
reject('axial-cut', lambda c: c.AddMetal('bad').AddLinPoly(
    points=square, norm_dir=2, elevation=-2, length=3.1), 'polygon end cuts a leaf')
reject('transformed', lambda c: c.AddMetal('bad').AddLinPoly(
    points=square, norm_dir=2, elevation=-2, length=4).AddTransform('Translate', [1, 0, 0]),
    'untransformed finite')
reject('noncartesian', lambda c: c.AddMetal('bad').AddLinPoly(
    points=square, norm_dir=2, elevation=-2, length=4).SetCoordinateSystem(1), 'Cartesian coordinates')
reject('material-polygon', lambda c: c.AddMaterial('bad', epsilon=2).AddLinPoly(
    points=square, norm_dir=2, elevation=-2, length=4), 'material objects must use axis-aligned volume boxes')
reject('rotated-body', lambda c: c.AddMetal('bad').AddRotPoly(
    points=square, norm_dir=2, elevation=0, rot_axis=1, angle=[0, np.pi]),
    'only planar or linearly extruded')
print('PASS: polygon-specific shape, coordinates, axial and material guards')

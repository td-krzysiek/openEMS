# SPDX-License-Identifier: GPL-3.0-or-later
"""Physical one-edge 50-ohm ports and local U/I integration.

Short microstrip fixture, identical copper and 50-ohm terminal coordinates.
The substrate cross section is intentionally coarse: this checks the adapter
and the refinement interface, not a calibrated 24-GHz manufacturing model.
"""
import json
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

Sim_Path = os.path.join(tempfile.gettempdir(), 'Cartesian_Subgrid_Ports')
os.makedirs(Sim_Path, exist_ok=True)
frequency = np.linspace(20e9, 28e9, 41)
results = {}


def fixture(kind, inside=False, reverse=False):
    options = dict(NrTS=24000, EndCriteria=0, TimeStepFactor=1, TimeStep=5e-14, OverSampling=20)
    if kind != 'yee':
        options['CartesianSubgridBox'] = [1.2, -.8, .2, 2.8, .8, 1.0]
        options['CartesianSubgridRatio'] = [1, 1, 1] if kind == 'identity' else ([2, 2, 1] if inside else [2, 1, 1])
    fdtd = openEMS(**options)
    fdtd.SetGaussExcite(24e9, 12e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(1e-3)
    fine = kind == 'yee'
    for axis, limits in zip('xyz', [(0, 4), (-1.2, 1.2), (0, 1.2)]):
        step = (.05 if fine and inside else .1) if axis != 'z' else .2
        lines = np.linspace(*limits, int(round((limits[1]-limits[0])/step))+1)
        if fine and not inside and axis == 'x':
            # Preserve the point-port cell, refine the propagation interval.
            # Globally shrinking a zero-width port also changes its local
            # capacitance/contact representation: that is a different fixture.
            lines = np.unique(np.round(np.r_[lines, np.arange(1.2, 2.8+.025, .05)], 12))
        mesh.SetLines(axis, np.round(lines, 12))
    metal = csx.AddMetal('electrodes')
    metal.AddBox([0, -1.2, .2], [4, 1.2, .2], priority=10)
    left, right = (1.4, 2.6) if inside else (.6, 3.4)
    metal.AddBox([left, -.2, .4], [right, .2, .4], priority=10)
    csx.AddMaterial('substrate', epsilon=2.25).AddBox([0, -1.2, .2], [4, 1.2, .4])
    ports = []
    for number, x in enumerate([left, right], 1):
        lo, hi = [x, 0, .2], [x, 0, .4]
        if reverse:
            lo, hi = hi, lo
        ports.append(fdtd.AddLumpedPort(number, 50, lo, hi, 'z', 1 if number == 1 else 0))
    return fdtd, csx, ports


def run(kind, inside=False, reverse=False):
    label = kind + ('-inside' if inside else '-outside') + ('-reversed' if reverse else '')
    fdtd, csx, ports = fixture(kind, inside, reverse)
    path = os.path.join(Sim_Path, label)
    os.makedirs(path, exist_ok=True)
    fdtd.Write2XML(os.path.join(path, 'input.xml'))
    old = os.getcwd()
    try:
        fdtd.Run(path, engine='basic', numThreads=1, exact_endcriteria=True)
    finally:
        os.chdir(old)
    for port in ports:
        port.CalcPort(path, frequency)
    p, q = ports
    values = dict(z=p.uf_tot/p.if_tot, s11=p.uf_ref/p.uf_inc, s21=q.uf_ref/p.uf_inc)
    for value in values.values():
        assert np.isfinite(value).all(), (label, 'nonfinite port result')
    assert np.max(np.abs(p.uf_inc)) > 1e-20, label + ': no excitation'
    power_error = float(np.max(np.abs(np.abs(values['s11'])**2+np.abs(values['s21'])**2-1)))
    assert power_error < 1e-3, (label, 'passive two-port power balance', power_error)
    values['power_error'] = power_error
    results[label] = values
    return values


# Identity ratio exercises the native source, resistor and probe path without
# an interface error. Compare to ordinary Yee on exactly the same base lines.
# 'yee' uses the matched tensor reference mesh; this first control uses an explicit
# replacement below so geometry and point-port contour sizes agree exactly.
fdtd, csx, ports = fixture('yee')
for axis, limits in zip('xy', [(0, 4), (-1.2, 1.2)]):
    csx.GetGrid().SetLines(axis, np.round(np.linspace(*limits, int(round((limits[1]-limits[0])/.1))+1), 12))
path = os.path.join(Sim_Path, 'yee-base')
os.makedirs(path, exist_ok=True)
fdtd.Write2XML(os.path.join(path, 'input.xml'))
old = os.getcwd()
try:
    fdtd.Run(path, engine='basic', numThreads=1, exact_endcriteria=True)
finally:
    os.chdir(old)
for port in ports:
    port.CalcPort(path, frequency)
p, q = ports
reference = dict(z=p.uf_tot/p.if_tot, s11=p.uf_ref/p.uf_inc, s21=q.uf_ref/p.uf_inc)
identity = run('identity')
metrics = {}
for name in ['port_ut_1', 'port_it_1', 'port_ut_2', 'port_it_2']:
    regular = np.loadtxt(os.path.join(path, name), comments='%')
    compact = np.loadtxt(os.path.join(Sim_Path, 'identity-outside', name), comments='%')
    np.testing.assert_allclose(regular[:, 0], compact[:, 0], rtol=1e-12, atol=0)
    peak = float(np.max(np.abs(regular[:, 1])))
    assert peak > 1e-12, (name, 'no excited waveform')
    relative = float(np.max(np.abs(regular[:, 1]-compact[:, 1]))/peak)
    assert relative < 2e-4, (name, 'absolute source/probe waveform mismatch', relative)
    metrics['identity-waveform-'+name] = relative
for key in ['z', 's11', 's21']:
    error = float(np.max(np.abs(identity[key]-reference[key])))
    scale = 50 if key == 'z' else 1
    metrics['identity-'+key] = error
    assert error/scale < 2e-4, (key, error)
print('native identical-grid 50-ohm port errors:', metrics, flush=True)

for inside in [False, True]:
    fine = run('yee', inside)
    local = run('local', inside)
    label = 'inside' if inside else 'outside'
    for key in ['z', 's11', 's21']:
        error = float(np.max(np.abs(local[key]-fine[key])))
        metrics[label+'-'+key] = error
        # Interface regression gates: 1% of 50 ohms, 0.005 complex S.
        # Same local cells at the point port are essential to this comparison.
        assert error < (.5 if key == 'z' else .005), (label, key, error)
    phase = float(np.max(np.abs(np.angle(local['s21']/fine['s21'], deg=True))))
    metrics[label+'-phase-deg'] = phase
    assert phase < .1, (label, 'transmission phase', phase)
    metrics[label+'-power-balance'] = local['power_error']
    print('local/fine', label, {k: v for k, v in metrics.items() if k.startswith(label)}, flush=True)

# Direction reversal must flip both U and I, leaving S and impedance unchanged.
reversed_result = run('local', True, True)
for key in ['z', 's11', 's21']:
    np.testing.assert_allclose(reversed_result[key], results['local-inside'][key], rtol=1e-7, atol=1e-7)

# Fail explicitly on a port with width, unresolved length or reactive element.
for name, start, stop, reactive in [
        ('distributed', [.6, -.2, .2], [.6, .2, .4], False),
        ('unresolved', [.6, 0, .2], [.6, 0, .3], False),
        ('overlapping', [.6, 0, .2], [.6, 0, .4], False),
        ('reactive', [.8, 0, .2], [.8, 0, .4], True)]:
    fdtd, csx, _ = fixture('local')
    if reactive:
        csx.AddLumpedElement('unsupported-C', ny=2, R=50, C=1e-12).AddBox(start, stop)
    else:
        fdtd.AddLumpedPort(3, 50, start, stop, 'z', 0)
    old = os.getcwd()
    try:
        try:
            fdtd.Run(os.path.join(Sim_Path, 'reject-'+name), engine='basic', setup_only=True)
        except (RuntimeError, ValueError) as error:
            assert 'CartesianSubgrid' in str(error)
        else:
            raise AssertionError(name + ': unsupported port accepted')
    finally:
        os.chdir(old)

with open(os.path.join(Sim_Path, 'metrics.json'), 'w') as output:
    json.dump(metrics, output, indent=2)
print('PASS: native 50-ohm resistor/source, complex Z/S11/S21, phase, reversal and guards')

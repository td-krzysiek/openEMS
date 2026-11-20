# SPDX-License-Identifier: GPL-3.0-or-later
"""Device port normalization against the same scoped double CPU microstrip.

Geometry copied explicitly from Cartesian_Subgrid_Ports.py; each candidate XML
must exactly match its current CPU control. This is not a calibrated PCB model.
Default tests the explicit reference backend. OPENEMS_SUBGRID_DEVICE=metal or
cuda requires a separately checked physical build, without availability skips.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from _subgrid_device_support import run_native_device

Sim_Path = os.path.join(tempfile.gettempdir(),'Cartesian_Subgrid_Device_Ports')
os.makedirs(Sim_Path,exist_ok=True)
frequency = np.linspace(20e9,28e9,41)
physical = os.environ.get('OPENEMS_SUBGRID_DEVICE','')
assert physical in ('','metal','cuda')


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


def run(label,engine,kind,inside,reverse):
    fdtd,csx,ports = fixture(kind,inside,reverse)
    # CSXCAD assigns random display colours; fix those without changing RF inputs.
    for prop in csx.GetAllProperties():
        prop.SetColor([96,128,160])
    path = Path(Sim_Path,label+'-'+engine)
    path.mkdir(exist_ok=False)
    fdtd.Write2XML(str(path/'input.xml'))
    previous = os.getcwd()
    try:
        result = run_native_device(fdtd, path, engine, physical, numThreads=1, exact_endcriteria=True)
        assert result in (None,0), (label,engine,result)
    finally:
        os.chdir(previous)
    for port in ports:
        port.CalcPort(str(path),frequency)
    first,second = ports
    values = dict(z=first.uf_tot/first.if_tot,s11=first.uf_ref/first.uf_inc,s21=second.uf_ref/first.uf_inc)
    assert all(np.isfinite(v).all() for v in values.values())
    balance = float(np.max(abs(abs(values['s11'])**2+abs(values['s21'])**2-1)))
    assert balance < 1e-3,(label,engine,'power',balance)
    np.savez(path/'results.npz',**values)
    return path,values,balance


def compare(a,b,va,vb):
    assert hashlib.sha256((a/'input.xml').read_bytes()).digest()==hashlib.sha256((b/'input.xml').read_bytes()).digest(), 'different native input'
    errors = {k:float(np.max(abs(va[k]-vb[k]))) for k in va}
    assert errors['z'] < .01 and errors['s11'] < 2e-4 and errors['s21'] < 2e-4,errors
    for name in ('port_ut_1','port_it_1','port_ut_2','port_it_2'):
        left,right = [np.loadtxt(p/name,comments='%') for p in (a,b)]
        assert left.shape==right.shape and np.isfinite(right).all()
        np.testing.assert_allclose(left[:,0],right[:,0],rtol=1e-12,atol=0)
        peak = float(np.max(abs(left[:,1])))
        assert peak > 1e-12,'unexcited port'
        error = float(np.max(abs(left[:,1]-right[:,1]))/peak)
        assert error < 2e-4,(name,error)
        errors[name] = error
    errors['phaseDeg'] = float(np.max(abs(np.angle(vb['s21']/va['s21'],deg=True))))
    assert errors['phaseDeg'] < .01,errors
    print(b.name,errors,flush=True)
    return errors


metrics = {}
for label,kind,inside,reverse in [('identity','identity',False,False),('outside','local',False,False),
                                    ('inside','local',True,False),('reverse','local',True,True)]:
    a,va,power = run(label,'basic',kind,inside,reverse)
    for engine in (['gpu-reference','gpu'] if physical else ['gpu-reference']):
        b,vb,devicePower = run(label,engine,kind,inside,reverse)
        metrics[label+'-'+engine] = dict(errors=compare(a,b,va,vb),cpuPowerError=power,devicePowerError=devicePower)
    Path(Sim_Path,'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print('PASS: native compact device ports, complex Z/S11/S21, phase, power, polarity and identical inputs',flush=True)

# SPDX-License-Identifier: GPL-3.0-or-later
"""Native fixed domains: matched XML, actual backends, U/I and TD/FD fields.

CPU CI uses two explicit reference workers. Physical qualification supplies
OPENEMS_SUBGRID_DEVICE and, optionally, OPENEMS_SUBGRID_DOMAIN_WORKERS. Two
logical Metal workers are not two physical devices. DFT remains on the host
at complete processing batches; no external collective-pause API is admitted.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile
import numpy as np
from _subgrid_device_support import run_native_device

HERE = Path(__file__).resolve().parent
ROOT = Path(tempfile.gettempdir()) / 'Cartesian_Subgrid_Domains'
ROOT.mkdir(exist_ok=True)
PHYSICAL = os.environ.get('OPENEMS_SUBGRID_DEVICE', '')
assert PHYSICAL in ('', 'metal', 'cuda')
WORKERS = os.environ.get('OPENEMS_SUBGRID_DOMAIN_WORKERS', 'auto,auto')
assert len(WORKERS.split(',')) == 2
metrics = {}


def definitions(name, marker):
    source = HERE / name
    scope = {'__file__': str(source), '__name__': 'native_domain_fixture'}
    exec(compile(source.read_text().split(marker, 1)[0], str(source), 'exec'), scope)
    return scope


fields = definitions('Cartesian_Subgrid_Device.py', '\nresults = {}')
ports = definitions('Cartesian_Subgrid_Device_Ports.py', '\nmetrics = {}')


def run(fdtd, label, engine, domains=False, axis='x', cuts='0'):
    path = ROOT / label
    path.mkdir(exist_ok=False)
    for prop in fdtd.GetCSX().GetAllProperties():
        prop.SetColor([96, 128, 160])
    fdtd.Write2XML(str(path / 'input.xml'))
    options = dict(numThreads=1, exact_endcriteria=True)
    if domains:
        options.update(subgrid_domain_axis=axis, subgrid_domain_cuts=cuts,
                       subgrid_domain_workers='reference,reference' if engine == 'gpu-reference' else WORKERS)
    previous = os.getcwd()
    try:
        run_native_device(fdtd, path, engine, PHYSICAL, **options)
    finally:
        os.chdir(previous)
    return path


def same_input(a, b):
    assert hashlib.sha256((a/'input.xml').read_bytes()).digest() == hashlib.sha256((b/'input.xml').read_bytes()).digest()


for loss in [False, True]:
    def fixture():
        fdtd, csx = fields['fixture'](steps=512, source_x=0)
        fdtd.SetBoundaryCond(['PML_4']*6)
        if loss:
            csx.AddMaterial('interior_loss', epsilon=1, kappa=.01, sigma=120).AddBox([-2]*3, [2]*3)
        return fdtd
    baseline = run(fixture(), f'fields-loss{loss}-single', 'gpu-reference')
    for axis in 'xyz':
        engines = ['gpu-reference'] + (['gpu'] if PHYSICAL else [])
        for engine in engines:
            candidate = run(fixture(), f'fields-loss{loss}-{axis}-{engine}', engine, True, axis)
            same_input(baseline, candidate)
            metrics[candidate.name] = fields['compare'](baseline, candidate)

frequency = ports['frequency']
for inside in [False, True]:
    data = {}
    for label, engine, domains in [('single', 'gpu-reference', False), ('domains', 'gpu-reference', True)] + ([('physical', 'gpu', True)] if PHYSICAL else []):
        fdtd, csx, terminal = ports['fixture']('local', inside)
        path = run(fdtd, f'ports-inside{inside}-{label}', engine, domains, cuts='2')
        for port in terminal:
            port.CalcPort(str(path), frequency)
        first, second = terminal
        values = dict(z=first.uf_tot/first.if_tot, s11=first.uf_ref/first.uf_inc,
                      s21=second.uf_ref/first.uf_inc)
        assert all(np.isfinite(value).all() for value in values.values())
        data[label] = path, values
        if label != 'single':
            a, va = data['single']
            metrics[path.name] = ports['compare'](a, path, va, values)

(ROOT/'metrics.json').write_text(json.dumps(dict(status='PASS_NATIVE_FIXED_DOMAIN_TD_FD_UI_CONTROLS',
    physical=PHYSICAL, workers=WORKERS if PHYSICAL else 'reference,reference', cases=metrics,
    limits=['Small native PML/loss fields and PEC two-port controls only.',
            'No representative antenna, production scale, remote workers or collective external control.']), indent=2)+'\n')
print('PASS: native fixed domains, matched XML, TD/FD fields, PML/loss and RF U/I ports')

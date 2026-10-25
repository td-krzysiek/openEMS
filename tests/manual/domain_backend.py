# SPDX-License-Identifier: GPL-3.0-or-later
"""Native multi-domain physics versus one device; bounded public fixtures.

Generate XML with matching *baseline* Python bindings, then execute the new
native CLI in child processes. Compare all complex E/H components and a probe
on the two-domain cut. No private antenna or production material fit is used.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
import re

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


def model(path, physics):
    fdtd = openEMS(NrTS=768, EndCriteria=0)
    fdtd.SetGaussExcite(3e9, 2e9)
    if physics == 'steady-state':
        fdtd.SetSinusExcite(3e9)
    fdtd.SetBoundaryCond(['MUR' if physics == 'mur' else 'PML_4']*6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis in 'xyz':
        grid.AddLine(axis, np.linspace(-30, 30, 25))
    csx.AddMaterial('loss', epsilon=2.1, kappa=.015).AddBox([-12.5]*3, [12.5]*3)
    csx.AddExcitation('source', exc_type=0, exc_val=[.3, .5, 1]).AddBox(
        [-2.5, -2.5, -2.5], [2.5, 2.5, 2.5], priority=20)
    csx.AddProbe('probe', p_type=2).AddPoint([0, 5, 5])
    if physics in ['pec', 'opaque']:
        metal = csx.AddMetal('physical_conductor')
        if physics == 'opaque':
            metal.AddAttribute('ExperimentalOpaqueSurfaceV1', '1;10,1e-12')
        metal.AddBox([-7.5, -7.5, -10], [7.5, 10, -5], priority=10)
    if physics == 'rlc':
        csx.AddLumpedElement('load', ny=2, caps=True, R=50, L=2e-9, C=2e-12).AddBox(
            [0, 7.5, 5], [0, 7.5, 10])
    for field, kind in [('E', 10), ('H', 11)]:
        csx.AddDump('field_'+field, dump_type=kind, dump_mode=1,
                    file_type=1, frequency=[3e9, 5e9]).AddBox([-20]*3, [20]*3)
    fdtd.Write2XML(str(path))


def fields(folder):
    result = {}
    for path in sorted(folder.glob('field_*.h5')):
        with h5py.File(path, 'r') as stream:
            for name, dataset in stream['FieldData/FD'].items():
                result[path.name+'/'+name] = dataset[()]
    assert result
    return result


def compare(candidate, control):
    a, b = fields(candidate), fields(control)
    assert a.keys() == b.keys()
    error = 0.
    for key in a:
        assert a[key].shape == b[key].shape
        assert np.isfinite(a[key]).all() and np.isfinite(b[key]).all()
        peak = float(np.max(abs(b[key])))
        assert peak > 0
        delta = float(np.max(abs(a[key]-b[key])))/peak
        assert delta < 1e-5, (key, delta)
        error = max(error, delta)
    x, y = [np.loadtxt(folder/'probe', comments='%', ndmin=2) for folder in [candidate, control]]
    assert x.shape == y.shape
    np.testing.assert_allclose(x[:, 0], y[:, 0], rtol=0, atol=0)
    delta = float(np.max(abs(x[:, 1:]-y[:, 1:])))/float(np.max(abs(y[:, 1:])))
    assert delta < 1e-5, delta
    return dict(max_complex_field_peak_relative_error=error, probe_peak_relative_error=delta)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True)
    parser.add_argument('--engine', default='gpu-reference')
    parser.add_argument('--output', required=True)
    parser.add_argument('--device-counts', default='2,3')
    args = parser.parse_args()
    root = Path(args.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    results = {}
    counts = [int(value) for value in args.device_counts.split(',')]
    for physics in ['vacuum', 'pec', 'rlc', 'opaque']:
        for count in [1, *counts]:
            folder = root/(physics+'-'+str(count))
            folder.mkdir()
            if count == 1:
                model(folder/'model.xml', physics)
            else:
                shutil.copy2(root/(physics+'-1/model.xml'), folder/'model.xml')
            command = [str(Path(args.binary).resolve()), 'model.xml', '--engine='+args.engine,
                       '--numThreads=2', '--exact-endcriteria', '-vv']
            if count > 1:
                command += ['--gpu-devices='+','.join(map(str, range(count)))]
            with (folder/'solver.log').open('w') as log:
                run = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=120)
            assert run.returncode == 0, (folder/'solver.log').read_text()[-3000:]
            text = (folder/'solver.log').read_text()
            assert 'Time for 768 iterations' in text
            if physics == 'pec':
                assert int(re.search(r'Number of PEC edges\s*:\s*(\d+)', text)[1]) > 0
            if physics == 'opaque':
                assert int(re.search(r'EXPERIMENTAL_OPAQUE_SURFACE_V1 active: (\d+)', text)[1]) > 0
            if physics == 'rlc':
                assert int(re.search(r'Active cells\s*:\s*(\d+)', text)[1]) > 0
            if args.engine == 'gpu':
                assert 'all extensions run on the device' in text and 'host fallback' not in text.lower()
            record = dict(command=command, xml_sha256=hashlib.sha256((folder/'model.xml').read_bytes()).hexdigest())
            if count > 1:
                record.update(compare(folder, root/(physics+'-1')))
            results[folder.name] = record
            (root/'results.json').write_text(json.dumps(results, indent=2)+'\n')
            print(folder.name, record, flush=True)
    if args.engine == 'gpu':
        # Mur is supported on a single device, but has no domain remapping.
        # A non-NULL generic HIP extension is insufficient evidence of support.
        for physics in ['mur', 'steady-state']:
            folder = root/('unsupported-'+physics)
            folder.mkdir()
            model(folder/'model.xml', physics)
            command = [str(Path(args.binary).resolve()), 'model.xml', '--engine=gpu',
                       '--gpu-devices=0,1', '--numThreads=2', '-vv']
            with (folder/'solver.log').open('w') as log:
                run = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=120)
            assert run.returncode != 0
            assert 'unsupported extension for native device-domain execution' in (folder/'solver.log').read_text()
            results['unsupported-'+physics] = dict(rejected=True, exit_code=run.returncode)
        (root/'results.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()

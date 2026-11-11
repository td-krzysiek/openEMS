# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare two native runtimes on the same public inputs, including a patch.

This checks an upgrade, not physical antenna accuracy or convergence. The
reference build supplies the same backend and thread count. A corrected CFL
calculation can change dt: record it and reject a change above 1e-8 relative.
Complex field/probe errors must stay below 1e-5 of each reference peak; the
patch's complex reflection error must stay below 1e-4 over 4--6 GHz.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import h5py
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


def fixture_module():
    source = Path(__file__).with_name('domain_backend.py')
    spec = importlib.util.spec_from_file_location('domain_fixture', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def patch(path):
    fdtd = openEMS(NrTS=4000, EndCriteria=0)
    fdtd.SetGaussExcite(5e9, 3e9)
    fdtd.SetBoundaryCond(['PML_4'] * 6)
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    grid = csx.GetGrid()
    grid.SetDeltaUnit(1e-3)
    for axis in 'xy':
        grid.AddLine(axis, np.linspace(-30, 30, 41))
    grid.AddLine('z', np.unique(np.r_[np.linspace(-21, 21, 29), -1.5]))
    csx.AddMaterial('substrate', epsilon=3.66).AddBox([-21, -21, -1.5], [21, 21, 0])
    metal = csx.AddMetal('metal')
    metal.AddBox([-21, -21, -1.5], [21, 21, -1.5], priority=10)
    metal.AddBox([-7.5, -9, 0], [7.5, 9, 0], priority=10)
    port = fdtd.AddLumpedPort(1, 50, [-3, -1.5, -1.5], [-3, 1.5, 0],
                              'z', excite=1, priority=5)
    fdtd.CreateNF2FFBox(start=[-22.5, -22.5, -13.5], stop=[22.5, 22.5, 13.5],
                       frequency=[4e9, 5e9, 6e9])
    fdtd.Write2XML(str(path))
    return fdtd, port  # keep the CSX and port owners alive while extracting


def arrays(folder):
    result = {}
    for path in sorted(folder.glob('*.h5')):
        with h5py.File(path, 'r') as stream:
            if 'FieldData/FD' in stream:
                for name, dataset in stream['FieldData/FD'].items():
                    result[path.name + '/' + name] = dataset[()]
    assert result, folder
    return result


def compare(candidate, reference):
    left, right = arrays(candidate), arrays(reference)
    assert left.keys() == right.keys()
    worst = 0.
    nonzero = 0
    for name in left:
        a, b = left[name], right[name]
        assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        peak = float(abs(b).max())
        error = float(abs(a - b).max())
        if peak:
            nonzero += 1
            error /= peak
        assert error <= 1e-5, (name, error)
        worst = max(worst, error)
    assert nonzero
    return dict(max_complex_field_peak_relative_error=worst, nonzero_arrays=nonzero)


def run(binary, folder, engine, devices=None):
    command = [str(binary), 'model.xml', '--engine=' + engine, '--numThreads=2',
               '--exact-endcriteria', '-vv']
    if devices:
        command.append('--gpu-devices=' + devices)
    started = time.monotonic()
    with (folder / 'solver.log').open('x') as log:
        proc = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT,
                              timeout=180)
    text = (folder / 'solver.log').read_text()
    assert proc.returncode == 0, text[-3000:]
    dt = float(re.search(r'Timestep \(s\)\s*:\s*([\d.eE+-]+)', text)[1])
    timing = re.search(r'Time for (\d+) iterations with ([\d.eE+-]+) cells\s*:\s*([\d.eE+-]+) sec', text)
    assert timing
    if engine == 'gpu':
        assert 'all extensions run on the device' in text and 'host fallback' not in text
    return dict(command=command, timestep_seconds=dt, steps=int(timing[1]),
                fdtd_seconds=float(timing[3]), wall_seconds=time.monotonic() - started,
                xml_sha256=hashlib.sha256((folder / 'model.xml').read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--baseline-binary', type=Path, required=True)
    parser.add_argument('--candidate-binary', type=Path, required=True)
    parser.add_argument('--engine', default='gpu-reference')
    parser.add_argument('--candidate-devices')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fixture = fixture_module()
    results = {}
    for physics in ['vacuum', 'pec', 'rlc', 'opaque', 'lorentz', 'lossy-metal', 'patch']:
        reference, candidate = [root / (physics + '-' + suffix) for suffix in ['old', 'new']]
        reference.mkdir()
        candidate.mkdir()
        owners = patch(reference / 'model.xml') if physics == 'patch' else None
        if owners is None:
            fixture.model(reference / 'model.xml', physics)
        shutil.copy2(reference / 'model.xml', candidate / 'model.xml')
        old = run(args.baseline_binary.resolve(), reference, args.engine)
        new = run(args.candidate_binary.resolve(), candidate, args.engine, args.candidate_devices)
        assert new['steps'] == old['steps']
        relative_dt = abs(new['timestep_seconds'] / old['timestep_seconds'] - 1)
        assert relative_dt < 1e-8, relative_dt
        result = dict(baseline=old, candidate=new, relative_timestep_change=relative_dt,
                      **compare(candidate, reference))
        if owners:
            frequency = np.linspace(4e9, 6e9, 201)
            port = owners[1]
            port.CalcPort(str(reference), frequency)
            old_gamma = port.uf_ref / port.uf_inc
            port.CalcPort(str(candidate), frequency)
            new_gamma = port.uf_ref / port.uf_inc
            assert np.isfinite(old_gamma).all() and np.isfinite(new_gamma).all()
            error = float(abs(new_gamma - old_gamma).max())
            assert error < 1e-4, error
            np.savez(root / 'patch-reflection.npz', frequency_hz=frequency,
                     baseline=old_gamma, candidate=new_gamma)
            result['max_complex_reflection_error'] = error
        results[physics] = result
        (root / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        print(physics, json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

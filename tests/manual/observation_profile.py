"""Compare observation monitoring off/on on a supplied public XML fixture.

Requires numpy/h5py. Each run has unique outputs; no production input is changed.
Native FDTD times exclude setup. Host-wall scopes can overlap and are not summed
into a claimed percentage of runtime. The caller preserves runtime/input identity.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import time

import h5py
import numpy as np


def compare(folder, reference):
    errors = []
    probes = []
    for path in reference.glob('*.h5'):
        with h5py.File(path) as before, h5py.File(folder/path.name) as after:
            if 'FieldData/FD' not in before:
                continue
            for key, value in before['FieldData/FD'].items():
                a, b = after['FieldData/FD'][key][:], value[:]
                assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
                scale = max(float(abs(b).max()), 1e-30)
                error = float(abs(a-b).max())/scale
                assert error <= 1e-6, (path.name, key, error)
                errors.append(error)
    for path in reference.iterdir():
        if path.name != 'probe' and not path.name.startswith(('port_ut_', 'port_it_')):
            continue
        a, b = [np.loadtxt(p, comments='%', ndmin=2) for p in [folder/path.name, path]]
        assert a.shape == b.shape and np.array_equal(a[:, 0], b[:, 0])
        error = float(abs(a[:, 1:]-b[:, 1:]).max())/max(float(abs(b[:, 1:]).max()), 1e-30)
        assert error <= 1e-6, (path.name, error)
        probes.append(error)
    assert errors and probes
    return dict(max_field_relative_error=max(errors), max_probe_relative_error=max(probes))


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--engine', default='basic')
    parser.add_argument('--solver-option', action='append', default=[])
    parser.add_argument('--repeats', type=int, default=2)
    args = parser.parse_args()
    assert args.repeats >= 1
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    result = dict(status='RUNNING', input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(), runs=[])
    reference = None
    for index in range(args.repeats*2):
        # Repeat ABBA blocks so both modes appear first/last.
        enabled = index % 4 in [1, 2]
        folder = root/(str(index)+('-on' if enabled else '-off'))
        folder.mkdir()
        shutil.copy2(args.input, folder/'model.xml')
        env = dict(os.environ, OPENEMS_OBSERVATION_PROFILE=str(folder/'observations.jsonl') if enabled else '')
        command = [str(args.binary.resolve()), 'model.xml', '--engine='+args.engine,
                   '--numThreads=2', '--exact-endcriteria', '-vv', *args.solver_option]
        start = time.monotonic()
        with (folder/'solver.log').open('x') as log:
            process = subprocess.run(command, cwd=folder, env=env, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=300)
        text = (folder/'solver.log').read_text()
        assert process.returncode == 0, text[-2500:]
        timing = re.search(r'Time for (\d+) iterations with ([\d.eE+-]+) cells\s*:\s*([\d.eE+-]+) sec', text)
        assert timing, text[-2500:]
        record = dict(enabled=enabled, steps=int(timing[1]), fdtd_seconds=float(timing[3]),
                      wall_seconds=time.monotonic()-start)
        if reference is None:
            reference = folder
        else:
            record.update(compare(folder, reference))
        profile = folder/'observations.jsonl'
        if enabled:
            rows = [json.loads(line) for line in profile.read_text().splitlines()]
            assert rows[0]['type'] == 'header' and rows[0]['device_waits_added'] is False
            scopes = [r for r in rows if r['type'] == 'scope']
            assert scopes and all(np.isfinite(r['host_wall_seconds']) and r['host_wall_seconds'] >= 0 for r in scopes)
            stages = {r['stage'] for r in scopes}
            assert {'dft_extract', 'dft_cpu_accumulate'} <= stages, stages
            samples = [r for r in scopes if r['stage'] == 'dft_cpu_accumulate']
            for owner in {r['owner'] for r in samples}:
                selected = [r for r in samples if r['owner'] == owner]
                assert len(selected) > 1
                steps = [r['step'] for r in selected]
                assert all(b-a == selected[0]['sample_interval'] for a,b in zip(steps, steps[1:]))
                assert all(r['logical_bytes'] > 0 and r['frequencies'] > 0 for r in selected)
            record['profile'] = dict(bytes=profile.stat().st_size, scopes=len(scopes),
                stages={stage: dict(calls=sum(r['stage']==stage for r in scopes),
                    host_wall_seconds=sum(r['host_wall_seconds'] for r in scopes if r['stage']==stage),
                    logical_bytes=sum(r['logical_bytes'] for r in scopes if r['stage']==stage)) for stage in sorted(stages)})
        else:
            assert not profile.exists()
        result['runs'].append(record)
        (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    off = [r['fdtd_seconds'] for r in result['runs'] if not r['enabled']]
    on = [r['fdtd_seconds'] for r in result['runs'] if r['enabled']]
    result.update(status='PASS', median_fdtd_off_seconds=statistics.median(off),
        median_fdtd_on_seconds=statistics.median(on),
        observed_overhead_fraction=statistics.median(on)/statistics.median(off)-1,
        caveat='Small-fixture timing noise; not antenna overhead or DMA/device timing.')
    (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='runs'}))


if __name__ == '__main__':
    main()

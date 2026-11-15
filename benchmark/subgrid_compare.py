#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Interleaved compact-core comparison; does not measure ordinary Yee/GPU speed.

Build subgrid_core.cpp against the previous and candidate core sources with
the same compiler/options. Candidate builds may define SUBGRID_FACTOR_STATS.
No Python solver bindings or extra packages are required by this driver.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(a, b):
    for name in ['leaves', 'edges', 'faces', 'steps', 'size', 'case']:
        assert a[name] == b[name], name
    assert abs(a['stable_timestep']/b['stable_timestep']-1) < 2e-13
    peak = max(abs(v) for v in a['fields'])
    assert peak > 1e-12 and len(a['fields']) == len(b['fields'])
    error = max(abs(x-y) for x, y in zip(a['fields'], b['fields']))/peak
    assert error < 2e-11, ('E/H dense reference mismatch', error)
    assert max(a['invariant_drift'], b['invariant_drift']) < 1e-10
    return error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert args.repeat >= 3, 'At least three repetitions are needed for the comparison'
    assert not args.output.exists(), 'Preserve prior measurements'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    binaries = {k: getattr(args, k).resolve() for k in ['baseline', 'candidate']}
    result = dict(utc=datetime.now(timezone.utc).isoformat(), host=platform.platform(),
        binaries={k: dict(path=str(p), sha256=sha(p)) for k, p in binaries.items()},
        scope='Standalone serial compact core, not native processing, ordinary Yee or GPU performance.',
        repeat=args.repeat, cases=[])
    cases = [('uniform', 16, 1000), ('local', 24, 500), ('multiple', 16, 500),
             ('strip', 64, 1000), ('strip', 128, 1000), ('strip', 256, 500), ('strip', 512, 200)]
    print('| Case | Size | Setup A/B [s] | Update A/B [s] | Accounted A/B [MiB] | Peak RSS A/B [MiB] | E/H error |', flush=True)
    print('|---|---:|---:|---:|---:|---:|---:|', flush=True)
    for kind, size, steps in cases:
        runs = {k: [] for k in binaries}
        errors = []
        for repeat in range(args.repeat):
            order = list(binaries) if repeat % 2 == 0 else list(reversed(binaries))
            current = {}
            for label in order:
                command = [str(binaries[label]), kind, str(size), str(steps)]
                completed = subprocess.run(command, capture_output=True, text=True, timeout=300)
                assert completed.returncode == 0, (command, completed.stderr)
                data = json.loads(completed.stdout)
                runs[label].append(data)
                current[label] = data
            errors.append(compare(current['baseline'], current['candidate']))
        metrics = {}
        for name in ['setup_seconds', 'update_seconds', 'accounted_bytes', 'peak_rss_bytes']:
            values = {label: [r[name] for r in rows] for label, rows in runs.items()}
            medians = {label: statistics.median(v) for label, v in values.items()}
            metrics[name] = dict(median=medians,
                spread_relative={k: (max(v)-min(v))/medians[k] if medians[k] else None for k, v in values.items()})
        row = dict(case=kind, size=size, steps=steps, runs=runs, metrics=metrics,
                   maximum_field_error=max(errors))
        result['cases'].append(row)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        a = lambda key: metrics[key]['median']['baseline']
        b = lambda key: metrics[key]['median']['candidate']
        print('| {} | {} | {:.4f}/{:.4f} | {:.4f}/{:.4f} | {:.2f}/{:.2f} | {:.2f}/{:.2f} | {:.2g} |'.format(
            kind, size, a('setup_seconds'), b('setup_seconds'), a('update_seconds'), b('update_seconds'),
            a('accounted_bytes')/2**20, b('accounted_bytes')/2**20,
            a('peak_rss_bytes')/2**20, b('peak_rss_bytes')/2**20, max(errors)), flush=True)
    result['status'] = 'PASS_FIELDS_TIMESTEP_INVARIANT_INTERLEAVED_MEASUREMENTS'
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()

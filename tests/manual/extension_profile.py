"""Bounded public aggregate/detail domain-profile comparison, with no extra waits.

Supply a public XML fixture and a matching installed core/plugin. Profiling
cadence is identical in aggregate and detail runs; only labels/output differ.
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

import h5py
import numpy as np


def compare(folder, reference):
    fields, probes = [], []
    for path in reference.glob('*.h5'):
        with h5py.File(path) as before, h5py.File(folder/path.name) as after:
            if 'FieldData/FD' not in before:
                continue
            for key, value in before['FieldData/FD'].items():
                a, b = after['FieldData/FD'][key][:], value[:]
                assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
                scale = float(abs(b).max())
                assert scale > 1e-30
                delta = float(abs(a-b).max())/scale
                assert delta <= 1e-6, (path.name, key, delta)
                fields.append(delta)
    for path in reference.iterdir():
        if path.name != 'probe' and not path.name.startswith(('port_ut_', 'port_it_')):
            continue
        a, b = [np.loadtxt(p, comments='%', ndmin=2) for p in [folder/path.name, path]]
        assert a.shape == b.shape and np.array_equal(a[:, 0], b[:, 0])
        scale = float(abs(b[:, 1:]).max())
        assert scale > 1e-30
        delta = float(abs(a[:, 1:]-b[:, 1:]).max())/scale
        assert delta <= 1e-6, (path.name, delta)
        probes.append(delta)
    assert fields and probes
    return dict(max_field_relative_error=max(fields), max_probe_relative_error=max(probes))


def check_profile(path, detailed):
    records = [json.loads(line) for line in path.read_text().splitlines()]
    assert records[0]['extension_detail'] is detailed
    assert records[0]['added_device_waits'] is False
    assert records[-1]['type'] == 'closed'
    samples = [r for r in records if r['type'] == 'sample']
    assert samples
    complete = 0
    for sample in samples:
        for domain in sample['domains']:
            if not detailed:
                assert 'extensions' not in domain
                continue
            rows = domain['extensions']
            assert rows and [r['id'] for r in rows] == list(range(len(rows)))
            for extension in rows:
                assert extension['name'] and len(extension['hooks']) == 7
                assert extension['calls'] == sum(h['calls'] for h in extension['hooks'])
                if extension['missing_device_intervals'] or not extension['calls']:
                    assert extension['device_seconds'] is None
                else:
                    assert extension['device_seconds'] >= 0
                    assert np.isclose(extension['device_seconds'],
                        sum(h['device_seconds'] or 0 for h in extension['hooks']), rtol=1e-6, atol=1e-12)
            if domain['device_timing_complete']:
                complete += 1
                assert np.isclose(sum(r['device_seconds'] or 0 for r in rows),
                                  domain['extension_device_seconds'], rtol=1e-6, atol=1e-12)
    if detailed and all('HIP' in d['backend'] for d in records[0]['domains']):
        assert complete > 0, 'HIP controls require usable device events'
    return dict(samples=len(samples), complete_detail_domain_samples=complete,
                interval_counts=[[d['interval_count'] for d in s['domains']] for s in samples])


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--devices', default='0,1')
    parser.add_argument('--cuts', required=True)
    parser.add_argument('--repeats', type=int, default=2)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    result = dict(status='RUNNING', source_xml_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(), runs=[])
    reference = None
    counts = None
    for index in range(-1, args.repeats*2):
        enabled = index % 4 in [1, 2] if index >= 0 else False
        folder = root/('plain' if index < 0 else str(index)+('-detail' if enabled else '-aggregate'))
        folder.mkdir()
        shutil.copy2(args.input, folder/'model.xml')
        env = dict(os.environ, OPENEMS_DOMAIN_PROFILE='' if index < 0 else str(folder/'profile.jsonl'),
                   OPENEMS_DOMAIN_PROFILE_EXTENSIONS='1' if enabled else '0',
                   OPENEMS_DOMAIN_PROFILE_WARMUP='128', OPENEMS_DOMAIN_PROFILE_INTERVAL='128')
        command = [str(args.binary.resolve()), 'model.xml', '--engine=gpu',
                   '--gpu-devices='+args.devices, '--gpu-domain-cuts='+args.cuts,
                   '--numThreads=2', '--exact-endcriteria', '-vv']
        with (folder/'solver.log').open('x') as log:
            process = subprocess.run(command, cwd=folder, env=env, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=300)
        text = (folder/'solver.log').read_text()
        assert process.returncode == 0, text[-2500:]
        timing = re.search(r'Time for (\d+) iterations with ([\d.eE+-]+) cells\s*:\s*([\d.eE+-]+) sec', text)
        assert timing
        record = dict(mode='plain' if index < 0 else ('detail' if enabled else 'aggregate'),
                      fdtd_seconds=float(timing[3]), steps=int(timing[1]))
        if reference is None:
            reference = folder
        else:
            record.update(compare(folder, reference))
        if index >= 0:
            record['profile'] = check_profile(folder/'profile.jsonl', enabled)
            if counts is None:
                counts = record['profile']['interval_counts']
            assert record['profile']['interval_counts'] == counts, 'labels added or removed event intervals'
        result['runs'].append(record)
        (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    aggregate = [r['fdtd_seconds'] for r in result['runs'] if r['mode']=='aggregate']
    detail = [r['fdtd_seconds'] for r in result['runs'] if r['mode']=='detail']
    result.update(status='PASS', aggregate_median_seconds=statistics.median(aggregate),
                  detail_median_seconds=statistics.median(detail),
                  observed_detail_overhead_fraction=statistics.median(detail)/statistics.median(aggregate)-1)
    (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='runs'}))


if __name__ == '__main__':
    main()

# SPDX-License-Identifier: GPL-3.0-or-later
"""Bounded native before/after checks; never modifies a production case."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import time
import xml.etree.ElementTree as ET


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def traces(folder):
    result = {}
    for name in ['p1_ub', 'p2_ub', 'shift-1_p1_ub', 'shift1_p1_ub']:
        rows = [list(map(float, line.split())) for line in (folder/name).read_text().splitlines()
                if line and not line.startswith('%')]
        assert rows and all(math.isfinite(v) for row in rows for v in row)
        result[name] = rows
    return result


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--base-bin', type=Path, required=True)
    p.add_argument('--patched-bin', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--backend', required=True)
    p.add_argument('--trace')
    p.add_argument('--large', action='store_true')
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    source = a.source.resolve()
    binaries = {'base': a.base_bin.resolve(), 'patched': a.patched_bin.resolve()}
    identity = {k: sha(v) for k, v in binaries.items()}
    result = {'backend': a.backend, 'binary_sha256': identity, 'runs': []}

    def run(label, variant, payload, decay=False):
        folder = a.output/label
        folder.mkdir()
        (folder/'model.xml').write_bytes(payload)
        env = dict(os.environ, OPENEMS_GPU_BACKEND=a.backend)
        if a.trace:
            env['LD_PRELOAD'] = a.trace
        start = time.monotonic()
        with (folder/'solver.log').open('w') as log:
            process = subprocess.run([str(binaries[variant]), 'model.xml', '--engine=gpu',
                                      '--numThreads=4', '--exact-endcriteria', '-vv'],
                                     cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600)
        text = (folder/'solver.log').read_text()
        assert process.returncode == 0, text[-3000:]
        assert 'all extensions run on the device' in text and 'backend: '+a.backend in text
        assert 'host fallback' not in text
        match = re.search(r'Time for (\d+) iterations with ([\deE+.-]+) cells : ([\deE+.-]+) sec', text)
        assert match, text[-3000:]
        if decay:
            assert 'end-criteria of -70.00dB reached' in text and int(match[1]) < 160000
        record = dict(run=label, variant=variant, steps=int(match[1]), locations=float(match[2]),
                      iteration_seconds=float(match[3]), wall_seconds=time.monotonic()-start,
                      XML_sha256=sha(folder/'model.xml'), log_sha256=sha(folder/'solver.log'))
        if a.trace:
            audit = re.search(r'CUDA_ALLOCATION_AUDIT calls=(\d+) peak_bytes=(\d+) remaining_bytes=(\d+)', text)
            assert audit and int(audit[3]) == 0, text[-3000:]
            record.update(cuda_calls=int(audit[1]), peak_bytes=int(audit[2]), remaining_bytes=int(audit[3]))
        data = traces(folder)
        result['runs'].append(record)
        (a.output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps(record), flush=True)
        return data, record

    payload = (source/'TESTSUITE/opaquesurface/graded_microstrip.xml').read_bytes()
    reference = None
    records = []
    # Interleaved ABBA controls simple drift; identical XML and exact energy schedule.
    for label, variant in [('a1','base'), ('b1','patched'), ('b2','patched'), ('a2','base')]:
        data, record = run('graded-'+label, variant, payload, True)
        assert all(max(abs(row[1]) for row in rows) > 1e-12 for rows in data.values())
        if reference is None:
            reference = data
        assert data == reference, 'Before/after full-decay waveform changed'
        records.append(record)
    assert len({r['steps'] for r in records}) == 1
    result['graded_iteration_speed_ratio'] = statistics.median(r['iteration_seconds'] for r in records if r['variant']=='base') / statistics.median(r['iteration_seconds'] for r in records if r['variant']=='patched')
    if a.trace:
        assert max(r['peak_bytes'] for r in records if r['variant']=='patched') < min(r['peak_bytes'] for r in records if r['variant']=='base')
    if a.large:
        root = ET.fromstring(payload)
        root.find('FDTD').set('NumberOfTimesteps','6000')
        root.find('FDTD').set('endCriteria','0')
        grid = root.find('.//RectilinearGrid')
        # Public physical line, refined grid only: not a qualified antenna model.
        for axis, factor in [('X',4),('Y',2)]:
            element = grid.find(axis+'Lines')
            values = list(map(float, element.text.split(',')))
            fine = [values[i]+(values[i+1]-values[i])*j/factor
                    for i in range(len(values)-1) for j in range(factor)] + [values[-1]]
            element.text = ','.join(format(v,'.17g') for v in fine)
            element.set('Qty',str(len(fine)))
        payload = ET.tostring(root, encoding='utf-8', xml_declaration=True)
        reference = None
        for label, variant in [('a1','base'), ('b1','patched'), ('b2','patched'), ('a2','base')]:
            data, record = run('large-'+label, variant, payload)
            if reference is None:
                reference = data
            assert data == reference, 'Before/after bounded large-grid waveform changed'
            assert record['steps']==6000
        large = [r for r in result['runs'] if r['run'].startswith('large-')]
        result['large_iteration_speed_ratio'] = statistics.median(r['iteration_seconds'] for r in large if r['variant']=='base') / statistics.median(r['iteration_seconds'] for r in large if r['variant']=='patched')
        result['large_scope'] = '6000-step public refined line; not energy convergence or full-antenna RF qualification'
    result['identical_native_probe_samples'] = True
    result['status'] = 'PASS'
    assert identity == {k:sha(v) for k,v in binaries.items()}
    (a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__':
    main()

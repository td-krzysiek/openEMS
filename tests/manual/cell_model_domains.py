"""Public field/probe controls for the CellModel/domain integration; no antenna inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
import h5py
import numpy as np


def data(folder):
    values = {}
    for file in sorted(folder.glob('field_*.h5')):
        with h5py.File(file) as f:
            for name, ds in f['FieldData/FD'].items():
                values[file.name + '/' + name] = ds[()]
    assert values
    return values


def compare(folder, reference):
    a, b = data(folder), data(reference)
    assert a.keys() == b.keys()
    errors = {}
    for key in a:
        assert a[key].shape == b[key].shape and np.isfinite(a[key]).all()
        peak = float(np.max(abs(b[key])))
        error = float(np.max(abs(a[key] - b[key])))/max(peak, 1e-30)
        assert error < 2e-5, (folder, key, error)
        errors[key] = error
    x, y = [np.loadtxt(p/'probe', comments='%', ndmin=2) for p in [folder, reference]]
    assert x.shape == y.shape and np.array_equal(x[:, 0], y[:, 0])
    error = float(np.max(abs(x[:, 1:] - y[:, 1:]))) / max(float(np.max(abs(y[:, 1:]))), 1e-30)
    assert error < 2e-5, (folder, 'probe', error)
    return dict(max_field_peak_relative_error=max(errors.values()), probe_peak_relative_error=error)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--binary', required=True)
    p.add_argument('--fixtures', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--engine', default='gpu-reference')
    p.add_argument('--physics', default='vacuum,pec,rlc,lorentz,lossy-metal')
    p.add_argument('--plugin')
    p.add_argument('--basic', action='store_true')
    p.add_argument('--fusion', default='off', choices=['off', 'apply', 'all'])
    p.add_argument('--halo', default='serial', choices=['serial', 'async', 'overlap'])
    p.add_argument('--two-only', action='store_true')
    p.add_argument('--default-tiling', action='store_true')
    p.add_argument('--library')
    a = p.parse_args()
    root = Path(a.output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(OPENEMS_GPU_TBLOCK='0', OPENEMS_GPU_FUSED='0', OPENEMS_DOMAIN_HALOS=a.halo, OPENEMS_DOMAIN_HOST_BUFFERS='copy' if a.halo == 'serial' else 'borrowed', OPENEMS_GPU_LOSSY_FUSION=a.fusion)
    if a.default_tiling:
        env.pop('OPENEMS_GPU_TBLOCK', None)
        env.pop('OPENEMS_GPU_FUSED', None)
    if a.plugin:
        env['OPENEMS_PLUGIN_PATH'] = str(Path(a.plugin).resolve())
    result = {'binary_sha256': hashlib.sha256(Path(a.binary).read_bytes()).hexdigest(), 'cases': {}}
    if a.library:
        result['library_sha256'] = hashlib.sha256(Path(a.library).read_bytes()).hexdigest()
    for physics in a.physics.split(','):
        cases = [('one', None), ('two', '0,12,25'), ('unequal', '0,7,25'), ('thin-three', '0,1,13,25')]
        if a.two_only:
            cases = cases[:3]
        if a.basic:
            cases = [('basic', None)] + cases
        for label, cuts in cases:
            folder = root/(physics+'-'+label)
            folder.mkdir()
            shutil.copyfile(Path(a.fixtures)/(physics+'.xml'), folder/'model.xml')
            engine = 'basic' if label == 'basic' else a.engine
            command = [str(Path(a.binary).resolve()), 'model.xml', '--engine='+engine, '--numThreads=2', '--exact-endcriteria', '-vv']
            if cuts:
                command.extend(['--gpu-devices='+','.join(map(str, range(len(cuts.split(','))-1))), '--gpu-domain-cuts='+cuts])
            started = time.monotonic()
            with (folder/'solver.log').open('w') as log:
                process = subprocess.run(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
            text = (folder/'solver.log').read_text()
            assert process.returncode == 0, (folder, text[-4000:])
            assert 'Time for 768 iterations' in text
            record = dict(command=command, wall_seconds=time.monotonic()-started, xml_sha256=hashlib.sha256((folder/'model.xml').read_bytes()).hexdigest())
            if cuts or (label == 'one' and a.basic):
                record.update(compare(folder, root/(physics+('-basic' if label == 'one' else '-one'))))
            if physics == 'lossy-metal':
                assert 'Operator_Ext_LossyMetal:' in text and 'fallback to PEC' not in text
            result['cases'][folder.name] = record
            (root/'results.json').write_text(json.dumps(result, indent=2)+'\n')
            print(folder.name, record, flush=True)

if __name__ == '__main__':
    main()

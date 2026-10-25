# SPDX-License-Identifier: GPL-3.0-or-later
"""Check one/two-device energy stops with the native cross-device error gate.

Partitioned and fused kernels need not be bit-identical. Use the same1e-5
reference-peak normalization as domain_backend.py; keep startup_gate.py's
stricter same-execution-path comparison unchanged.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np

LIMIT = 1e-5


def peak_relative_error(candidate, reference):
    assert candidate.shape == reference.shape
    assert np.isfinite(candidate).all() and np.isfinite(reference).all()
    peak = float(np.max(abs(reference)))
    difference = float(np.max(abs(candidate - reference)))
    if peak == 0:
        assert difference == 0, 'Nonzero result for a zero reference'
        return 0.0
    error = difference / peak
    assert error < LIMIT, error
    return error


def compare(candidate, control, fixture):
    a, b = fixture.fields(candidate), fixture.fields(control)
    assert a.keys() == b.keys()
    errors = [peak_relative_error(a[key], b[key]) for key in a]
    assert any(np.max(abs(values)) > 0 for values in b.values())
    x, y = [np.loadtxt(folder / 'probe', comments='%', ndmin=2)
            for folder in [candidate, control]]
    assert x.shape == y.shape
    np.testing.assert_array_equal(x[:, 0], y[:, 0])
    return dict(max_complex_field_peak_relative_error=max(errors),
                nonzero_arrays=sum(int(np.max(abs(values)) > 0) for values in b.values()),
                probe_peak_relative_error=peak_relative_error(x[:, 1:], y[:, 1:]),
                reference_peak_relative_limit=LIMIT,
                comparison_basis='Same native cross-device gate as domain_backend.py; exact probe times')


def main():
    spec = importlib.util.spec_from_file_location('fixture', Path(__file__).with_name('startup_gate.py'))
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True)
    parser.add_argument('--two-device-control', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    root, control = Path(args.output).resolve(), Path(args.two_device_control).resolve()
    root.mkdir(parents=True, exist_ok=False)
    shutil.copy2(control / 'model.xml', root / 'model.xml')
    with (root / 'solver.log').open('w') as log:
        subprocess.run([str(Path(args.binary).resolve()), 'model.xml', '--engine=gpu',
                        '--gpu-device=0', '--numThreads=2', '--exact-endcriteria', '-v'],
                       cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=120)
    pattern = r'end-criteria of -50(?:\.0+)?dB reached after (\d+) timesteps'
    one = int(re.search(pattern, (root / 'solver.log').read_text())[1])
    two = int(re.search(pattern, (control / 'solver.log').read_text())[1])
    assert one == two, (one, two)
    result = dict(one_gpu_steps=one, two_gpu_steps=two, **compare(root, control, fixture))
    (root / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

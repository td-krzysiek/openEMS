# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 openEMS contributors
"""Public synthetic timing fixture: 201x201 Hertzian-dipole recording planes.

No FDTD or private antenna inputs. Reports read + transform + write wall time
at 3 GHz on a 61x120 angular grid, four CPU threads versus strict GPU. Requires
a build with NF2FF_HIP and an available device. Run in an isolated GPU slot.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from openEMS import _nf2ff
from openEMS.nf2ff import nf2ff_results
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'python/Tests'))
import test_nf2ff_backend as fixture


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert _nf2ff._nf2ff.GetGpuDevice(), 'A GPU build and device are required'
    fixture.N_LINES = 201
    fixture._write_box(str(out))
    theta = np.deg2rad(np.arange(0, 181, 3))
    phi = np.deg2rad(np.arange(0, 360, 3))
    values, timings = {}, []
    for label, backend in [('cpu', 'cpu'), ('gpu-1', 'gpu'), ('gpu-2', 'gpu')]:
        start = time.monotonic()
        transform = _nf2ff._nf2ff([fixture.FREQ], theta, phi, [0, 0, 0], numThreads=4)
        transform.SetBackend(backend)
        for face in range(6):
            assert transform.AnalyseFile(str(out / ('nf2ff_E_%d.h5' % face)),
                                         str(out / ('nf2ff_H_%d.h5' % face)))
        transform.Write2HDF5(str(out / (label + '.h5')))
        timings.append(dict(label=label, seconds=time.monotonic() - start))
        values[label] = nf2ff_results(str(out / (label + '.h5')))
    cpu = values['cpu']
    peak = float(np.max(np.abs(cpu.E_norm[0])))
    checks = {}
    for label in ['gpu-1', 'gpu-2']:
        gpu = values[label]
        error = max(float(np.max(np.abs(a[0] - b[0]))) / peak for a, b in
                    [(gpu.E_theta, cpu.E_theta), (gpu.E_phi, cpu.E_phi)])
        assert error < 1e-4
        np.testing.assert_allclose(gpu.Dmax, cpu.Dmax, rtol=1e-4)
        np.testing.assert_allclose(gpu.Prad, cpu.Prad, rtol=1e-6)
        checks[label] = dict(relative_complex_field_error=error)
    result = dict(device=_nf2ff._nf2ff.GetGpuDevice(), frequencies_hz=[fixture.FREQ],
                  surface_lines=201, angular_shape=[len(theta), len(phi)],
                  cpu_threads=4, timings=timings, checks=checks)
    (out / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

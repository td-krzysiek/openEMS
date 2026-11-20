# SPDX-License-Identifier: GPL-3.0-or-later
"""Shared native log/backend checks for scoped Cartesian-subgrid device tests."""
import ctypes
import os
import re
from pathlib import Path
import sys


def run_native_device(fdtd, path, engine, physical='', **options):
    """Capture native diagnostics and reject backend substitution or setup failure."""
    if physical not in ('', 'metal', 'cuda'):
        raise ValueError('unknown physical subgrid test selector')
    if engine == 'gpu' and not physical:
        raise ValueError('an explicit physical qualification selector is required')
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    if (path / 'solver.log').exists():
        raise FileExistsError('preserve prior native evidence; use a new case directory')
    options.setdefault('verbose', 1)
    if engine == 'gpu' and physical == 'cuda':
        options['gpu_device'] = 0
    library = ctypes.CDLL('msvcrt' if os.name == 'nt' else None)
    library.fflush.argtypes = [ctypes.c_void_p]
    library.fflush.restype = ctypes.c_int
    sys.stdout.flush()
    sys.stderr.flush()
    library.fflush(None)
    saved = [os.dup(1), os.dup(2)]
    try:
        with (path / 'solver.log').open('w') as log:
            os.dup2(log.fileno(), 1)
            os.dup2(log.fileno(), 2)
            try:
                result = fdtd.Run(str(path), engine=engine, **options)
            finally:
                sys.stdout.flush()
                sys.stderr.flush()
                library.fflush(None)
                os.dup2(saved[0], 1)
                os.dup2(saved[1], 2)
    finally:
        for descriptor in saved:
            os.close(descriptor)
    text = (path / 'solver.log').read_text(errors='replace')
    assert result in (None, 0), (path.name, result, text[-2000:])
    if engine in ('gpu', 'gpu-reference'):
        expected = 'reference (' if engine == 'gpu-reference' else ('Metal (' if physical == 'metal' else 'HIP (')
        if 'subgrid_domain_workers' in options:
            workers = options['subgrid_domain_workers'].split(',')
            selected = re.findall(r'^CartesianSubgrid native domain (\d+): (.+)$', text, re.MULTILINE)
            assert len(selected) == len(workers), ('missing native workers', path.name, selected)
            for index, ((owner, backend), worker) in enumerate(zip(selected, workers)):
                required = 'reference (' if worker == 'reference' else expected
                assert int(owner) == index and backend.startswith(required), (path.name, owner, backend, worker)
            assert 'host mirror at complete processing batches' in text
        else:
            assert 'CartesianSubgrid device: ' + expected in text, ('wrong native device backend', path.name, text[-2000:])
        if engine == 'gpu' and physical == 'cuda':
            assert re.search(r'^\s*0: .+ \(hip-nvidia, ', text, re.MULTILINE), (
                'CUDA requested but selected device0 is not the NVIDIA HIP plugin', path.name, text[-2000:])
    return result

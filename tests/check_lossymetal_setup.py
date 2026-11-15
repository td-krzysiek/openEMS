#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Full setup-state regression. Standard library only; no Python bindings.

python3 tests/check_lossymetal_setup.py build/tests/inspect_lossymetal_setup /tmp/lossy-check
Optional third argument: an inspector linked against the unchanged parent.
The output directory must be new. This test never starts an FDTD engine.
"""
from pathlib import Path
import hashlib
import os
import shutil
import subprocess
import sys

helper = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
reference = Path(sys.argv[3]).resolve() if len(sys.argv) > 3 else helper
out.mkdir(parents=True, exist_ok=False)
fixtures = Path(__file__).parent / 'fixtures/lossymetal-setup'
for source in sorted(fixtures.glob('*.xml')):
    expected = None
    variants = [(reference, 'serial', 4)]
    variants += [(helper, mode, threads) for mode, threads in
                 [('serial', 4), ('cache', 4), ('parallel-no-cache', 4),
                  ('parallel', 1), ('parallel', 2), ('parallel', 4)]]
    for index, (binary, mode, threads) in enumerate(variants):
        folder = out / source.stem / f'{index}-{mode}-{threads}'
        folder.mkdir(parents=True)
        shutil.copy2(source, folder / 'model.xml')
        env = dict(os.environ, OPENEMS_LOSSY_SETUP=mode,
                   OPENEMS_PROFILE_LOSSY_SETUP='0')
        with (folder / 'setup.log').open('w') as log:
            result = subprocess.run([str(binary), 'model.xml', str(threads)] + ([] if source.stem in ('unchanged-pec', 'zero-height-no-thickness') else ['--expect-surface']),
                                    cwd=folder, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=180)
        if result.returncode:
            raise RuntimeError(f'Failed setup: {folder}')
        files = [folder / 'ade.bin', *folder.glob('operator_dump.*'),
                 *folder.glob('PEC_dump.*')]
        actual = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
        if not any(name.startswith('operator_dump') for name in actual):
            raise RuntimeError('Missing base coefficient dump')
        if expected is None:
            expected = actual
        if actual != expected:
            raise RuntimeError(f'ADE/base/PEC mismatch: {folder}')
    print('PASS', source.stem, flush=True)

# A deliberately disabled surface must be rejected, without reading an extension
# after CalcECOperator has deleted it. This reproduces the old empty-state/UAF.
folder = out / 'negative-excluded-surface'
folder.mkdir()
source = fixtures / 'separate-boundary-bodies.xml'
shutil.copy2(source, folder / 'model.xml')
with (folder / 'setup.log').open('w') as log:
    result = subprocess.run([str(helper), 'model.xml', '4', '--expect-surface', '--exclude-all'],
                            cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=60)
if result.returncode != 8:
    raise RuntimeError('Excluded surface was not rejected safely')
print('PASS excluded-surface negative control', flush=True)

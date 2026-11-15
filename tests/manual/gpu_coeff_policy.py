# SPDX-License-Identifier: GPL-3.0-or-later
"""Fixed-step public auto/raw coefficient equivalence, independently of fusion."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess

from gpu_lossy_fusion import model


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--logical-metal', action='store_true')
    p.add_argument('--layouts', default='single;0,9,25;0,12,13,25')
    args = p.parse_args()
    spec = importlib.util.spec_from_file_location('compare', Path(__file__).with_name('extension_profile.py'))
    comp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comp)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    rows = []
    for case in ['smooth', 'rough', 'sheet', 'curved-void', 'multiple', 'vacuum', 'lorentz']:
        for layout in args.layouts.split(';'):
            reference = None
            for name, policy, fusion, binary in [('parent', 'auto', 'off', args.baseline),
                    ('auto', 'auto', 'off', args.binary), ('raw', 'raw', 'off', args.binary),
                    ('raw-fused', 'raw', 'all', args.binary)]:
                folder = root/(case+'-'+layout.replace(',', '_')+'-'+name)
                folder.mkdir()
                model(folder/'model.xml', case)
                binary = binary.resolve()
                prefix = binary.parents[1]
                env = dict(os.environ, OPENEMS_GPU_MAIN_COEFF=policy, OPENEMS_GPU_LOSSY_FUSION=fusion,
                    OPENEMS_DOMAIN_HALOS='serial' if args.logical_metal else 'overlap',
                    OPENEMS_DOMAIN_HOST_BUFFERS='copy' if args.logical_metal else 'borrowed',
                    OPENEMS_PLUGIN_PATH=str(prefix.parent/'build-metal/test-plugins') if args.logical_metal
                        else str(prefix/'lib/openEMS/plugins'),
                    LD_LIBRARY_PATH=str(prefix/'lib')+':/reference/deps/lib:/opt/openEMS/lib')
                command = [str(binary), 'model.xml', '--engine=gpu', '--numThreads=2', '--exact-endcriteria', '-vv']
                if layout != 'single':
                    count = len(layout.split(','))-1
                    command += ['--gpu-devices='+','.join(str(i) for i in range(count)), '--gpu-domain-cuts='+layout]
                with (folder/'solver.log').open('x') as log:
                    r = subprocess.run(command, env=env, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=240)
                text = (folder/'solver.log').read_text()
                assert r.returncode == 0 and 'Time for 2048 iterations' in text, text[-3000:]
                assert 'all extensions run on the device' in text
                assert ('GPU main coefficient policy: raw' in text) == (policy == 'raw'), folder
                enabled = fusion == 'all' and case not in ['vacuum', 'lorentz']
                assert ('GPU LossyMetal branch fusion: all' in text) == enabled, folder
                row = dict(case=case, layout=layout, variant=name)
                if reference is None: reference = folder
                else: row.update(comp.compare(folder, reference))
                rows.append(row)
                print(folder.name, 'PASS', flush=True)
                (root/'result.json').write_text(json.dumps(dict(status='RUNNING', runs=rows), indent=2)+'\n')
    (root/'result.json').write_text(json.dumps(dict(status='PASS', tolerance=1e-6, runs=rows), indent=2)+'\n')


if __name__ == '__main__':
    main()

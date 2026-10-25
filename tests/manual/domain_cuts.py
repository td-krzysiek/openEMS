# SPDX-License-Identifier: GPL-3.0-or-later
"""Native uneven X ownership versus frozen unsplit public fixtures.

Include one-node edge/interior slabs: their halo-to-owned ratio and cut PML/
source/material locations catch errors hidden by equal slabs. The Metal test
plugin exercises real kernels on one physical GPU; this is not a speed test.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location('fixture', Path(__file__).with_name('domain_backend.py'))
fixture = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixture)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--binary', required=True)
    parser.add_argument('--engine', default='gpu-reference')
    parser.add_argument('--control', required=True, help='Completed domain_backend.py output')
    parser.add_argument('--output', required=True)
    parser.add_argument('--max-domains', type=int, choices=[2, 3], default=3)
    args = parser.parse_args()
    binary, control, root = map(lambda p: Path(p).resolve(), [args.binary, args.control, args.output])
    root.mkdir(parents=True, exist_ok=False)
    results = {}
    plans = [[0, 1, 25], [0, 24, 25], [0, 8, 25]]
    if args.max_domains == 3:
        plans += [[0, 1, 24, 25], [0, 8, 9, 25]]
    for physics in ['vacuum', 'pec', 'rlc', 'opaque']:
        for cuts in plans:
            folder = root/(physics+'-'+'-'.join(map(str, cuts)))
            folder.mkdir()
            shutil.copy2(control/(physics+'-1')/'model.xml', folder/'model.xml')
            grid = ET.parse(folder/'model.xml').find('ContinuousStructure/RectilinearGrid/XLines')
            assert grid is not None and len(grid.text.split(',')) == 25
            command = [str(binary), 'model.xml', '--engine='+args.engine, '--numThreads=2', '-vv',
                       '--exact-endcriteria', '--gpu-devices='+','.join(map(str, range(len(cuts)-1))),
                       '--gpu-domain-cuts='+','.join(map(str, cuts))]
            with (folder/'solver.log').open('w') as log:
                run = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=120)
            text = (folder/'solver.log').read_text()
            assert run.returncode == 0 and 'Time for 768 iterations' in text, text[-3000:]
            for n in range(len(cuts)-1):
                assert f'Domain {n}: owns X[{cuts[n]},{cuts[n+1]})' in text
            if args.engine == 'gpu':
                assert 'all extensions run on the device' in text and 'host fallback' not in text.lower()
            results[folder.name] = dict(xml_sha256=hashlib.sha256((folder/'model.xml').read_bytes()).hexdigest(),
                                        command=command, **fixture.compare(folder, control/(physics+'-1')))
            (root/'results.json').write_text(json.dumps(results, indent=2)+'\n')
            print(folder.name, results[folder.name], flush=True)
    # Fail closed on malformed or incomplete ownership; never quietly use equal cuts.
    invalid = ['0,25', '0,0,25', '0,15,14', '1,12,25', '0,12,24', '0,12,26',
               '0,-1,25', '0,1.5,25', '0,,25', '0,12,25,', '', '0,4294967296,25']
    for n, cuts in enumerate(invalid + [None]):
        folder = root/('invalid-'+str(n))
        folder.mkdir()
        shutil.copy2(control/'vacuum-1/model.xml', folder/'model.xml')
        command = [str(binary), 'model.xml', '--engine='+args.engine, '--numThreads=1', '-vv']
        if cuts is not None:
            command += ['--gpu-devices=0,1', '--gpu-domain-cuts='+cuts]
        else:
            command += ['--gpu-domain-cuts=0,12,25']
        with (folder/'solver.log').open('w') as log:
            run = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=30)
        text = (folder/'solver.log').read_text()
        assert run.returncode != 0 and 'Running FDTD engine' not in text, text[-3000:]
        results[folder.name] = dict(rejected=True, cuts=cuts, exit_code=run.returncode)
    (root/'results.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()

# SPDX-License-Identifier: GPL-3.0-or-later
"""Cell-model GPU domain equivalence on independent public XML fixtures.

Runs one executor and unequal slabs, including one-line slabs in a conductor,
source and both physical boundaries. Compares all complex E/H components and
probe samples; the reference backend exercises the actual coordinator on CPU.
Use --engine=gpu with a matching physical GPU plugin for hardware validation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

from gpu_lossy_fusion import model as conductor_model
from extension_profile import compare


def model(path, case, steps):
    conductor_model(path, case if case not in ['opaque', 'rlc'] else 'smooth', steps)
    tree = ET.parse(path)
    tree.getroot().find('FDTD').set('endCriteria', '0')
    props = tree.getroot().find('.//Properties')
    if case == 'opaque':
        metal = props.find('LossyMetal')
        metal.tag = 'Metal'
        metal.attrib = dict(ID='3', Name='opaque')
        ET.SubElement(metal, 'Attributes', ExperimentalOpaqueSurfaceV1='1;10,1e-12')
    elif case == 'rlc':
        props.remove(props.find('LossyMetal'))
        load = ET.SubElement(props, 'LumpedElement', ID='3', Name='load', Direction='2',
                             Caps='1', R='50', C='2e-12', L='2e-9', LEtype='0')
        box = ET.SubElement(ET.SubElement(load, 'Primitives'), 'Box', Priority='0')
        ET.SubElement(box, 'P1', X='0', Y='7.5', Z='5')
        ET.SubElement(box, 'P2', X='0', Y='7.5', Z='10')
    tree.write(path, encoding='utf-8', xml_declaration=True)


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--engine', default='gpu-reference')
    p.add_argument('--cases', default='vacuum,smooth,rough,sheet,curved-void,multiple,lorentz,opaque,rlc')
    p.add_argument('--layouts', default='single;0,9,25;0,12,13,25;0,1,24,25')
    p.add_argument('--steps', type=int, default=768)
    p.add_argument('--halos', choices=['serial', 'async', 'overlap'], default='serial')
    args = p.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve()
    env = dict(os.environ, OPENEMS_DOMAIN_HALOS=args.halos)
    for name in ['OPENEMS_GPU_LOSSY_FUSION', 'OPENEMS_GPU_MAIN_COEFF']:
        if env.get(name):
            p.error('unset legacy test control '+name)
    records = []
    for case in args.cases.split(','):
        reference = None
        for layout in args.layouts.split(';'):
            folder = root/(case+'-'+layout.replace(',', '_'))
            folder.mkdir()
            model(folder/'model.xml', case, args.steps)
            cmd = [str(binary), 'model.xml', '--engine='+args.engine,
                   '--numThreads=2', '--exact-endcriteria', '-vv']
            if layout != 'single':
                count = len(layout.split(','))-1
                cmd += ['--gpu-devices='+','.join(map(str, range(count))), '--gpu-domain-cuts='+layout]
            with (folder/'solver.log').open('x') as log:
                run = subprocess.run(cmd, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
            text = (folder/'solver.log').read_text()
            assert run.returncode == 0, (case, layout, run.returncode, text[-2500:])
            assert f'Time for {args.steps} iterations' in text
            if case == 'opaque':
                assert 'EXPERIMENTAL_OPAQUE_SURFACE_V1 active:' in text
            elif case not in ['vacuum', 'lorentz', 'rlc']:
                assert 'Surface edges' in text and 'fallback to PEC' not in text
            record = dict(case=case, layout=layout, command=cmd,
                          xml_sha256=hashlib.sha256((folder/'model.xml').read_bytes()).hexdigest())
            if reference is None:
                reference = folder
            else:
                record.update(compare(folder, reference))
            records.append(record)
            (root/'result.json').write_text(json.dumps(dict(status='RUNNING', runs=records), indent=2)+'\n')
            print(folder.name, 'PASS', flush=True)
    (root/'result.json').write_text(json.dumps(dict(status='PASS', tolerance=1e-6, runs=records), indent=2)+'\n')


if __name__ == '__main__':
    main()

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
    conductor_model(path, 'sheet' if case=='positive-sheet' else
                    (case if case not in ['opaque', 'rlc'] else 'smooth'), steps)
    tree = ET.parse(path)
    tree.getroot().find('FDTD').set('endCriteria', '0')
    props = tree.getroot().find('.//Properties')
    if case == 'positive-sheet':
        metal = props.find('LossyMetal')
        metal.tag = 'ConductingSheet'
        metal.attrib = dict(ID='3', Name='positive-sheet', Conductivity='1e4', Thickness='1e-5')
        ET.SubElement(metal, 'Attributes', ExperimentalSheetAdmittanceV1='1;10,1e-12;20,2e-12')
    elif case == 'opaque':
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



def check_legacy_controls(binary, root, environment):
    checks = []
    for name, value, message in [
            ('OPENEMS_GPU_MAIN_COEFF', 'raw', 'only supports auto'),
            ('OPENEMS_GPU_LOSSY_FUSION', 'all', 'unset OPENEMS_GPU_LOSSY_FUSION')]:
        folder = root/('rejected-'+name.lower())
        folder.mkdir()
        model(folder/'model.xml', 'vacuum', 32)
        with (folder/'solver.log').open('x') as log:
            run = subprocess.run([str(binary), 'model.xml', '--engine=gpu-reference', '--numThreads=2'],
                                 cwd=folder, env=dict(environment, **{name: value}), stdout=log,
                                 stderr=subprocess.STDOUT, timeout=30)
        assert run.returncode != 0 and message in (folder/'solver.log').read_text()
        checks.append(dict(control=name, value=value, rejected=True))
    return checks


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--engine', default='gpu-reference')
    p.add_argument('--cases', default='vacuum,smooth,rough,sheet,curved-void,multiple,lorentz,opaque,rlc,positive-sheet')
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
            elif case == 'positive-sheet':
                assert 'EXPERIMENTAL_POSITIVE_RL_SHEET_V1 active:' in text
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
    controls = check_legacy_controls(binary, root, env)
    (root/'result.json').write_text(json.dumps(dict(status='PASS', tolerance=1e-6, runs=records,
                                                   rejected_legacy_controls=controls), indent=2)+'\n')


if __name__ == '__main__':
    main()

# SPDX-License-Identifier: GPL-3.0-or-later
"""Fixed-step public GPU LossyMetal fusion equivalence, without CSXCAD bindings.

Use matching installed core/backends. Defaults cover a single GPU and unequal
two-/three-domain cuts, including a one-owned-line domain through the metal.
Fields and probe traces must agree within 1e-6 of each control's nonzero peak.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def model(path, physics, steps=2048):
    fixture = 'lorentz-control.xml' if physics == 'lorentz' else 'lossymetal-branches.xml'
    tree = ET.parse(Path(__file__).parent/'fixtures'/fixture)
    tree.getroot().find('FDTD').set('NumberOfTimesteps', str(steps))
    # Native XML treats a zero criterion as "use the default", not "disabled".
    tree.getroot().find('FDTD').set('endCriteria', '1e-80')
    props = tree.getroot().find('.//Properties')
    metal = props.find('LossyMetal')
    if physics == 'lorentz':
        pass
    elif physics == 'vacuum':
        props.remove(metal)
    elif physics == 'rough':
        metal.set('RoughnessModel', 'gradient')
        metal.set('Rq', '4e-7')
    elif physics == 'sheet':
        metal.set('Conductivity', '1e4')
        metal.set('Thickness', '1e-5')
        box = metal.find('.//Box')
        box.find('P2').set('Z', box.find('P1').get('Z'))
    elif physics == 'curved-void':
        primitives = metal.find('Primitives')
        primitives.clear()
        cylinder = ET.SubElement(primitives, 'Cylinder', Priority='10', Radius='9')
        ET.SubElement(cylinder, 'P1', X='0', Y='0', Z='-10')
        ET.SubElement(cylinder, 'P2', X='0', Y='0', Z='-5')
        air = ET.SubElement(props, 'Material', ID='6', Name='air-drill', Isotropy='1')
        ET.SubElement(air, 'Property', Epsilon='1', Mue='1')
        holes = ET.SubElement(air, 'Primitives')
        hole = ET.SubElement(holes, 'Cylinder', Priority='30', Radius='2.5')
        ET.SubElement(hole, 'P1', X='0', Y='0', Z='-12.5')
        ET.SubElement(hole, 'P2', X='0', Y='0', Z='-2.5')
    elif physics == 'multiple':
        extra = copy.deepcopy(metal)
        extra.set('ID', '6')
        extra.set('Name', 'second-conductivity')
        extra.set('Conductivity', '2e6')
        extra.find('.//P1').set('Z', '7.5')
        extra.find('.//P2').set('Z', '10')
        props.append(extra)
    elif physics != 'smooth':
        raise ValueError(physics)
    tree.write(path, encoding='utf-8', xml_declaration=True)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--layouts', default='single;0,9,25;0,12,13,25')
    parser.add_argument('--cases', default='smooth,rough,sheet,curved-void,multiple,vacuum,lorentz')
    parser.add_argument('--engine', default='gpu')
    parser.add_argument('--logical-metal', action='store_true', help='Use each build prefix sibling build-metal/test-plugins')
    parser.add_argument('--steps', type=int, default=2048)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('comparison', Path(__file__).with_name('extension_profile.py'))
    comp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(comp)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    runs = []
    for case in args.cases.split(','):
        for layout in args.layouts.split(';'):
            reference = None
            variants = [('off', args.binary), ('apply', args.binary), ('all', args.binary)]
            if args.baseline:
                variants.insert(0, ('parent', args.baseline))
            for variant, binary in variants:
                folder = out/(case+'-'+layout.replace(',', '_')+'-'+variant)
                folder.mkdir()
                model(folder/'model.xml', case, args.steps)
                binary = binary.resolve()
                prefix = binary.parents[1]
                env = dict(os.environ, OPENEMS_GPU_LOSSY_FUSION='off' if variant == 'parent' else variant,
                           OPENEMS_DOMAIN_HALOS='overlap', OPENEMS_DOMAIN_HOST_BUFFERS='borrowed',
                           LD_LIBRARY_PATH=str(prefix/'lib')+':/reference/deps/lib:/opt/openEMS/lib',
                           OPENEMS_PLUGIN_PATH=str(prefix/'lib/openEMS/plugins'))
                if args.logical_metal:
                    plugins = prefix.parent/'build-metal/test-plugins'
                    assert plugins.is_dir(), plugins
                    env['OPENEMS_PLUGIN_PATH'] = str(plugins)
                if args.logical_metal or args.engine == 'gpu-reference':
                    env['OPENEMS_DOMAIN_HALOS'] = 'serial'
                    env['OPENEMS_DOMAIN_HOST_BUFFERS'] = 'copy'
                command = [str(binary), 'model.xml', '--engine='+args.engine, '--numThreads=2',
                           '--exact-endcriteria', '-vv']
                if layout != 'single':
                    count = len(layout.split(','))-1
                    command += ['--gpu-devices='+','.join(str(i) for i in range(count)),
                                '--gpu-domain-cuts='+layout]
                with (folder/'solver.log').open('x') as log:
                    result = subprocess.run(command, cwd=folder, env=env, stdout=log,
                                            stderr=subprocess.STDOUT, timeout=240)
                text = (folder/'solver.log').read_text()
                assert result.returncode == 0, text[-3000:]
                assert f'Time for {args.steps} iterations' in text, folder
                if args.engine == 'gpu':
                    assert 'all extensions run on the device' in text, folder
                    enabled = variant in ['apply', 'all'] and case not in ['vacuum', 'lorentz']
                    assert ('GPU LossyMetal branch fusion: '+variant in text) == enabled, folder
                if case not in ['vacuum', 'lorentz']:
                    assert 'Surface edges' in text and 'fallback to PEC' not in text, folder
                record = dict(case=case, layout=layout, variant=variant,
                              xml_sha256=hashlib.sha256((folder/'model.xml').read_bytes()).hexdigest())
                if reference is None:
                    reference = folder
                else:
                    record.update(comp.compare(folder, reference))
                runs.append(record)
                (out/'result.json').write_text(json.dumps(dict(status='RUNNING', runs=runs), indent=2)+'\n')
                print(folder.name, 'PASS', flush=True)
    (out/'result.json').write_text(json.dumps(dict(status='PASS', tolerance=1e-6, runs=runs), indent=2)+'\n')


if __name__ == '__main__':
    main()

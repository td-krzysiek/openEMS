# SPDX-License-Identifier: GPL-3.0-or-later
"""Interleaved fixed-step public copper/storage benchmark; no production inputs.

This synthetic graded grid is not an antenna speedup prediction. Profiled hook
attribution must be run separately from these uninstrumented timing samples.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import time
import xml.etree.ElementTree as ET

from gpu_lossy_fusion import model


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--shape', default='193,257,97')
    p.add_argument('--steps', type=int, default=1536)
    p.add_argument('--repeats', type=int, default=4)
    p.add_argument('--devices', default='0')
    p.add_argument('--cuts', help='optional fixed X cuts, applied identically to both builds')
    p.add_argument('--fusion-only', action='store_true')
    p.add_argument('--baseline-fusion', choices=['off', 'apply', 'all'], default='off')
    p.add_argument('--one-step', action='store_true', help='disable temporal blocking, whole-step fusion and automatic trials on both builds')
    args = p.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    shape = [int(n) for n in args.shape.split(',')]
    assert len(shape) == 3 and min(shape) >= 25 and args.steps > 0 and args.repeats >= 2
    model(root/'model.xml', 'smooth', args.steps)
    tree = ET.parse(root/'model.xml')
    grid = tree.getroot().find('.//RectilinearGrid')
    for axis,n,power in zip('XYZ', shape, [1.031, 1.057, 1.0]):
        node = grid.find(axis+'Lines')
        node.set('Qty', str(n))
        node.text = ','.join(format(-30+60*(i/(n-1))**power, '.17g') for i in range(n))
    props = tree.getroot().find('.//Properties')
    metal = props.find('LossyMetal')
    box = metal.find('.//Box')
    for point in [box.find('P1'), box.find('P2')]:
        point.set('X', '-24' if point.tag == 'P1' else '24')
        point.set('Y', '-24' if point.tag == 'P1' else '24')
    extra = copy.deepcopy(metal)
    extra.set('ID', '6'); extra.set('Name', 'upper-copper')
    extra.find('.//P1').set('Z', '7.5'); extra.find('.//P2').set('Z', '10')
    props.append(extra)
    # Keep the full-update/metal/PML workload; bounded surface dumps retain an
    # independent complex-field check without making CPU volume DFT the target.
    for dump in props.findall('DumpBox'):
        dump.find('.//P1').set('Z', '0'); dump.find('.//P2').set('Z', '0')
    tree.write(root/'model.xml', encoding='utf-8', xml_declaration=True)
    raw = (root/'model.xml').read_bytes()
    spec = importlib.util.spec_from_file_location('comparison', Path(__file__).with_name('extension_profile.py'))
    comp = importlib.util.module_from_spec(spec); spec.loader.exec_module(comp)
    variants = [('parent', args.baseline, args.baseline_fusion, 'auto'), ('off', args.binary, 'off', 'auto'),
                ('apply', args.binary, 'apply', 'auto'), ('all', args.binary, 'all', 'auto')]
    if not args.fusion_only:
        variants += [('raw', args.binary, 'off', 'raw'), ('raw-all', args.binary, 'all', 'raw')]
    result = dict(status='RUNNING', synthetic=True, input_sha256=hashlib.sha256(raw).hexdigest(),
                  shape=shape, steps=args.steps, devices=args.devices, runs=[],
                  baseline_fusion=args.baseline_fusion, one_step=args.one_step, cuts=args.cuts)
    reference = None
    for repeat in range(args.repeats):
        # Alternate forward/reverse order to reduce systematic thermal/order bias.
        for name,binary,fusion,policy in variants if repeat%2 == 0 else reversed(variants):
            folder = root/(f'{repeat:02d}-'+name); folder.mkdir()
            (folder/'model.xml').write_bytes(raw)
            binary = binary.resolve(); prefix = binary.parents[1]
            env = {k:v for k,v in os.environ.items() if k not in ['OPENEMS_DOMAIN_PROFILE',
                'OPENEMS_DOMAIN_PROFILE_EXTENSIONS', 'OPENEMS_HIP_KERNEL_TIMES', 'OPENEMS_OBSERVATION_PROFILE']}
            env.update(OPENEMS_GPU_LOSSY_FUSION=fusion, OPENEMS_GPU_MAIN_COEFF=policy,
                OPENEMS_DOMAIN_HALOS='overlap', OPENEMS_DOMAIN_HOST_BUFFERS='borrowed',
                OPENEMS_PLUGIN_PATH=str(prefix/'lib/openEMS/plugins'),
                LD_LIBRARY_PATH=str(prefix/'lib')+':/reference/deps/lib:/opt/openEMS/lib')
            if args.one_step:
                env.update(OPENEMS_GPU_TBLOCK='0', OPENEMS_GPU_FUSED='0', OPENEMS_GPU_TRIAL='0')
            command = [str(binary), 'model.xml', '--engine=gpu', '--numThreads=4', '--exact-endcriteria', '-vv']
            if ',' in args.devices:
                command += ['--gpu-devices='+args.devices]
            else:
                command += ['--gpu-device='+args.devices]
            if args.cuts:
                assert ',' in args.devices, 'cuts require multiple devices'
                command += ['--gpu-domain-cuts='+args.cuts]
            started = time.monotonic()
            with (folder/'solver.log').open('x') as log:
                run = subprocess.run(command, env=env, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=900)
            text = (folder/'solver.log').read_text()
            assert run.returncode == 0, text[-3000:]
            match = re.search(r'Time for (\d+) iterations with ([\d.eE+-]+) cells : ([\d.eE+-]+) sec', text)
            assert match and int(match[1]) == args.steps, text[-2000:]
            assert 'all extensions run on the device' in text
            assert ('GPU LossyMetal branch fusion: '+fusion in text) == (fusion != 'off')
            assert ('GPU main coefficient policy: raw' in text) == (policy == 'raw')
            surface = re.search(r'Surface edges\s*:\s*(\d+)', text)
            assert surface and int(surface[1]) > 1000
            record = dict(variant=name, repeat=repeat, fdtd_seconds=float(match[3]),
                          process_seconds=time.monotonic()-started, surface_edges=int(surface[1]))
            if reference is None: reference = folder
            else: record.update(comp.compare(folder, reference))
            result['runs'].append(record)
            (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
            print(name, record['fdtd_seconds'], flush=True)
    medians = {name: statistics.median(r['fdtd_seconds'] for r in result['runs'] if r['variant']==name)
               for name,*_ in variants}
    result.update(status='PASS', median_fdtd_seconds=medians,
                  speedup_vs_parent={name: medians['parent']/value for name,value in medians.items()},
                  limitations=['Public synthetic workload; not full-antenna RF or speed qualification.',
                               'Median native FDTD includes observation/energy processing; setup is separate.'])
    (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()

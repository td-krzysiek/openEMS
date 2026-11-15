# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare copied/borrowed pinned halo payloads on public HIP inputs.

Requires installed CSXCAD/openEMS Python bindings for fixture generation. Use
physical GPUs exclusively; two logical devices on one GPU are correctness-only.
This never starts a project antenna or replaces a running runtime.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import time
import xml.etree.ElementTree as ET

import numpy as np
from domain_backend import model as fixture
from runtime_upgrade import patch, compare


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--devices', default='0,1')
    parser.add_argument('--benchmark', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve(strict=True)
    result = dict(status='RUNNING', binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                  devices=args.devices, controls={}, benchmarks=[])

    def save():
        (root/'results.json').write_text(json.dumps(result, indent=2)+'\n')

    def run(label, source, mode, cuts, host_mode='copy'):
        folder = root/label
        folder.mkdir()
        shutil.copy2(source, folder/'model.xml')
        env = dict(os.environ, OPENEMS_DOMAIN_HALOS=mode,
                   OPENEMS_DOMAIN_HOST_BUFFERS=host_mode,
                   OPENEMS_DOMAIN_PROFILE=str(folder/'profile.jsonl'),
                   OPENEMS_DOMAIN_PROFILE_WARMUP='128', OPENEMS_DOMAIN_PROFILE_INTERVAL='128')
        command = [str(binary), 'model.xml', '--engine=gpu', '--gpu-devices='+args.devices,
                   '--gpu-domain-cuts='+','.join(map(str,cuts)), '--numThreads=2',
                   '--exact-endcriteria', '-vv']
        started = time.monotonic()
        with (folder/'solver.log').open('x') as log:
            process = subprocess.run(command, cwd=folder, env=env, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=240)
        text = (folder/'solver.log').read_text()
        assert process.returncode == 0, (label,text[-2500:])
        assert 'all extensions run on the device' in text and 'host fallback' not in text.lower()
        assert 'Domain halo mode: '+mode in text and 'HIP (' in text
        assert 'Domain host buffers: '+host_mode in text
        assert 'Unused primitive' not in text and 'fallback to PEC' not in text
        timing = re.search(r'Time for (\d+) iterations with ([\d.eE+-]+) cells\s*:\s*([\d.eE+-]+) sec', text)
        assert timing, text[-2000:]
        records = [json.loads(line) for line in (folder/'profile.jsonl').read_text().splitlines()]
        samples = [r for r in records if r['type']=='sample']
        assert records[-1]['type']=='closed' and samples
        assert records[0]['halo_mode']==mode
        assert records[0]['host_buffer_mode']==host_mode
        assert all(d['device_timing_complete'] for s in samples for d in s['domains'])
        assert all(sum(d['halo_sent_bytes']for d in s['domains'])==sum(d['halo_received_bytes']for d in s['domains']) for s in samples)
        record = dict(mode=mode,host_mode=host_mode,cuts=cuts,steps=int(timing[1]),fdtd_seconds=float(timing[3]),
                      wall_seconds=time.monotonic()-started,
                      input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      median_exchange_scope_ms=1000*statistics.median(s['exchange_wall_seconds']for s in samples))
        time.sleep(.5)
        return folder, record

    inputs = root/'inputs'
    inputs.mkdir()
    for physics in ['vacuum','pec','rlc','opaque','lorentz','lossy-metal','patch']:
        source=inputs/(physics+'.xml')
        owners=patch(source) if physics=='patch' else fixture(source,physics)
        nx = int(ET.parse(source).getroot().find('ContinuousStructure/RectilinearGrid/XLines').get('Qty'))
        entries={}
        control=None
        for label,mode,host,cut in [('async-copy','async','copy',nx//2),
                                     ('async-borrowed','async','borrowed',nx//2),
                                     ('overlap-copy','overlap','copy',nx//2),
                                     ('overlap-borrowed','overlap','borrowed',nx//2),
                                     ('thin-borrowed','overlap','borrowed',1)]:
            folder,record=run(physics+'-'+label,source,mode,[0,cut,nx],host)
            if control is None:
                control=folder
            else:
                record.update(compare(folder,control))
                probe_errors=[]
                for p in control.iterdir():
                    if p.name=='probe' or p.name.startswith(('port_ut','port_it')):
                        a=np.loadtxt(folder/p.name,comments='%',ndmin=2)
                        b=np.loadtxt(p,comments='%',ndmin=2)
                        assert a.shape==b.shape and np.isfinite(a).all()
                        np.testing.assert_allclose(a[:,0],b[:,0],rtol=1e-8,atol=1e-20)
                        peak=float(abs(b[:,1:]).max())
                        assert peak>0
                        error=float(abs(a[:,1:]-b[:,1:]).max())/peak
                        assert error<=1e-5,(physics,label,p.name,error)
                        probe_errors.append(error)
                assert probe_errors
                record['max_probe_peak_relative_error']=max(probe_errors)
                if physics=='patch':
                    port=owners[1]
                    freq=np.linspace(4e9,6e9,201)
                    port.CalcPort(str(control),freq)
                    old=port.uf_ref/port.uf_inc
                    port.CalcPort(str(folder),freq)
                    new=port.uf_ref/port.uf_inc
                    assert np.isfinite(new).all()
                    error=float(abs(new-old).max())
                    assert error<1e-4,error
                    record['max_complex_reflection_error']=error
            entries[label]=record
            result['controls'][physics]=entries
            save()
            print(physics,label,record,flush=True)
    result['controls_status']='PASS'
    if args.benchmark:
        tree=ET.parse(inputs/'vacuum.xml')
        node=tree.getroot()
        node.find('FDTD').set('NumberOfTimesteps','1024')
        node.find('FDTD').set('endCriteria','0')
        for axis,count in zip('XYZ',[129,1535,106]):
            lines=node.find('ContinuousStructure/RectilinearGrid/'+axis+'Lines')
            lines.set('Qty',str(count))
            lines.text=','.join(format(-30+60*i/(count-1),'.17g')for i in range(count))
        props=node.find('ContinuousStructure/Properties')
        for child in list(props):
            if child.tag in ['DumpBox','ProbeBox']:
                props.remove(child)
        source=inputs/'benchmark.xml'
        tree.write(source)
        for trial in range(4):
            order = ['copy', 'borrowed'] if trial % 2 == 0 else ['borrowed', 'copy']
            for host in order:
                _,record=run('bench-'+str(trial)+'-'+host,source,'overlap',[0,65,129],host)
                result['benchmarks'].append(record)
                save()
        medians = {host: statistics.median(r['fdtd_seconds'] for r in result['benchmarks']
                   if r['host_mode'] == host) for host in ['copy', 'borrowed']}
        result['benchmark_median_fdtd_seconds'] = medians
        result['benchmark_speedup'] = medians['copy']/medians['borrowed']
        result['benchmark_scope'] = 'Same binary, input, partition and physical device pair; synthetic vacuum/PML, not antenna throughput'
    result['status']='PASS'
    save()


if __name__=='__main__':
    main()

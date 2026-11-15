#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Small numerical engine comparisons, using the native CLI and stdlib only.

OPENEMS_BIN=/path/to/openEMS OPENEMS_TEST_GPU=1 python3 test_gpu.py -v
Without OPENEMS_TEST_GPU, the basic/GPU-reference comparison still runs.
These coarse synthetic lines check indexing/state/order, not RF accuracy.
The separate test_stability.py checks decay on the original graded fixture.
"""
import cmath
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent


def fixture(rotation,pec=False,uniform=False):
    root=ET.parse(HERE/'graded_microstrip.xml').getroot()
    fdtd=root.find('FDTD')
    fdtd.set('NumberOfTimesteps','6000');fdtd.set('endCriteria','0')
    for face in fdtd.find('BoundaryCond').attrib:fdtd.find('BoundaryCond').set(face,'MUR')
    grid=root.find('ContinuousStructure/RectilinearGrid')
    lines=[[-2+i*.5 for i in range(49)],[-3+i*.5 for i in range(13)],
           [-1.5,-1.,-.5,-.035,0.,.13589,.27178,.30678,.5,1.,1.5,2.]]
    lines[1]+= [-.30060727534,0.,.30060727534]
    for node in root.iter():
        if all(a in node.attrib for a in 'XYZ'):
            for n,axis in enumerate('XYZ'):
                value=round(float(node.get(axis)),12)
                node.set(axis,str(value));lines[n].append(value)
    for n,axis in enumerate('XYZ'):
        element=grid.find(axis+'Lines')
        values=sorted(set(round(v,12) for v in lines[n]));element.text=','.join(map(str,values))
        element.set('Qty',str(len(values)))
    metals=root.findall('ContinuousStructure/Properties/Metal')
    if pec:
        for metal in metals:metal.remove(metal.find('Attributes'))
    else:
        # Different branch counts per fit; the third rotation exercises 128.
        metals[0].find('Attributes').set('ExperimentalOpaqueSurfaceV1','0.3;0.5,2e-13')
        if rotation==2:
            fit='0.3'+''.join(';{},{}'.format(.2+.001*j,(j+1)*1e-13) for j in range(128))
            metals[1].find('Attributes').set('ExperimentalOpaqueSurfaceV1',fit)
        if uniform:
            # Exercise the tiled device layout, including 128-pole histories.
            fit=metals[1].find('Attributes').get('ExperimentalOpaqueSurfaceV1')
            metals[0].find('Attributes').set('ExperimentalOpaqueSurfaceV1',fit)
    if rotation:
        # Permute all coordinates and component vectors together, including
        # material/excitation boxes, probes and mesh. All PML faces are equal.
        names=list('XYZ');axes=list(grid)
        for node in root.iter():
            if all(a in node.attrib for a in names):
                old=[node.get(a) for a in names]
                for n,a in enumerate(names):node.set(a,old[(n-rotation)%3])
            for attr in ['Excite','PropDir']:
                if attr in node.attrib:
                    old=node.get(attr).split(',')
                    node.set(attr,','.join(old[(n-rotation)%3] for n in range(3)))
        for n,node in enumerate(axes):node.tag=names[(n+rotation)%3]+'Lines'
    return root


def run(binary,folder,root,engine):
    folder.mkdir(parents=True)
    ET.ElementTree(root).write(folder/'model.xml')
    command=[binary,'model.xml','--engine='+engine,'--numThreads=4','--exact-endcriteria','-vv']
    with (folder/'solver.log').open('w') as log:
        result=subprocess.run(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=120)
    text=(folder/'solver.log').read_text()
    assert result.returncode==0,text[-3000:]
    assert 'Time for 6000 iterations' in text,text[-3000:]
    assert 'nan' not in text.lower() and ' inf' not in text.lower()
    if engine=='gpu':
        assert 'Engine_GPU: all extensions run on the device' in text,text[-3000:]
        assert 'host fallback' not in text,text[-3000:]
        backend=os.environ.get('OPENEMS_GPU_BACKEND')
        if backend:assert 'backend: '+backend in text
    traces={}
    for name in ['p1_ub','p2_ub','shift-1_p1_ub','shift1_p1_ub']:
        traces[name]=[list(map(float,line.split())) for line in (folder/name).read_text().splitlines()
                      if line and not line.startswith('%')]
        assert traces[name] and all(math.isfinite(x) for row in traces[name] for x in row)
        assert max(abs(row[1]) for row in traces[name])>1e-12,name+' is zero'
    return traces,text


@unittest.skipUnless(os.environ.get('OPENEMS_BIN'),'Set OPENEMS_BIN for compiled regressions')
class OpaqueGPU(unittest.TestCase):
    def test_rotated_boxes_variable_orders_and_pec(self):
        binary=os.path.abspath(os.environ['OPENEMS_BIN'])
        engines=['basic','gpu-reference']
        if os.environ.get('OPENEMS_TEST_GPU'):engines.append('gpu')
        with tempfile.TemporaryDirectory(prefix='openems-opaque-gpu-') as temp:
            output=os.environ.get('OPENEMS_TEST_OUTPUT')
            parent=Path(output)/'compact' if output else Path(temp)
            cases=[(0,False,False),(1,False,False),(2,False,False),(0,True,False)]
            cases += [(r,False,True) for r in range(3)]
            for rotation,pec,uniform in cases:
                root=fixture(rotation,pec,uniform);reference=None
                for engine in engines:
                    kind='pec' if pec else ('uniform' if uniform else 'surface')
                    folder=parent/('r{}-{}-{}'.format(rotation,kind,engine))
                    traces,text=run(binary,folder,root,engine)
                    if not pec:self.assertIn('EXPERIMENTAL_OPAQUE_SURFACE_V1 active:',text)
                    errors={}
                    if reference is None:reference=traces
                    for name,rows in traces.items():
                        ref=reference[name]
                        self.assertEqual(len(rows),len(ref))
                        self.assertEqual([r[0] for r in rows],[r[0] for r in ref])
                        peak=max(abs(r[1]) for r in ref)
                        error=max(abs(a[1]-b[1]) for a,b in zip(ref,rows))/peak
                        self.assertLess(error,1e-4,(engine,rotation,name,error))
                        complex_error=0.
                        for f in [22e9,24.125e9,26e9]:
                            a=sum(v*cmath.exp(-2j*math.pi*f*t) for t,v in ref)
                            b=sum(v*cmath.exp(-2j*math.pi*f*t) for t,v in rows)
                            delta=abs(a-b)/max(abs(a),1e-30)
                            self.assertLess(delta,5e-4,(engine,rotation,name,f,delta))
                            complex_error=max(complex_error,delta)
                        errors[name]={'peak_relative':error,'complex_relative':complex_error}
                    (folder/'comparison.json').write_text(json.dumps(dict(engine=engine,
                        rotation=rotation,pec=pec,uniform=uniform,errors=errors,
                        binary_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest()),indent=2)+'\n')


if __name__=='__main__':unittest.main()

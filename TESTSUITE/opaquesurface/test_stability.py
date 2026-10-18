#!/usr/bin/env python3
"""Opt-in compiled regression; standard-library Python only.

OPENEMS_BIN=/path/to/openEMS python3 test_stability.py
Use OPENEMS_BASELINE_BIN for an optional negative control: the unfixed
experimental boundary must trigger the energy-growth watchdog.
"""
import cmath
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import unittest
import xml.etree.ElementTree as ET

HERE=Path(__file__).resolve().parent


def solve(binary,engine='multithreaded',pec=False):
    with tempfile.TemporaryDirectory(prefix='openems-surface-cfl-') as temp:
        folder=Path(temp);root=ET.parse(HERE/'graded_microstrip.xml').getroot()
        if pec:
            for metal in root.findall('ContinuousStructure/Properties/Metal'):
                metal.remove(metal.find('Attributes'))
        ET.ElementTree(root).write(folder/'model.xml')
        env=os.environ.copy();env.update(OPENEMS_TILE_X='10',OPENEMS_TILE_Y='10',OPENEMS_TILE_STEPS='3')
        start=time.monotonic();reason=None
        with (folder/'solver.log').open('w') as log:
            process=subprocess.Popen([binary,'model.xml','--engine='+engine,'--numThreads=4',
                '--exact-endcriteria'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
            try:
                while process.poll() is None:
                    time.sleep(.5);text=(folder/'solver.log').read_text()
                    energies=[float(v) for v in re.findall(r'Energy: ~([\deE+.-]+)',text)]
                    if max(energies,default=0)>1e-5:reason='energy growth'
                    if time.monotonic()-start>1200:reason='timeout'
                    if reason:process.terminate();break
                process.wait(timeout=10)
            finally:
                if process.poll() is None:process.kill();process.wait()
        text=(folder/'solver.log').read_text()
        step=re.search(r'FDTD timestep is: ([\deE+.-]+)',text)
        samples=[];dfts={}
        for name in ['p1_ub','p2_ub']:
            path=folder/name
            if path.exists():
                samples.extend(float(v) for line in path.read_text().splitlines()
                    if line and not line.startswith('%') for v in line.split())
        for name in ['shift-1_p1_ub','p1_ub','shift1_p1_ub']:
            path=folder/name
            if path.exists():
                values=[list(map(float,line.split())) for line in path.read_text().splitlines()
                    if line and not line.startswith('%')]
                if len(values)>1:
                    dt=values[1][0]-values[0][0]
                    dfts[name]=sum(v*cmath.exp(-2j*math.pi*24.125e9*t) for t,v in values)*dt
        result=dict(exit=process.returncode,reason=reason,converged='reached after' in text,
            capped='Max. number of timesteps was reached' in text,
            dt=float(step.group(1)) if step else None,finite=all(math.isfinite(v) for v in samples),
            nonzero=any(abs(v)>1e-12 for v in samples[1::2]),dfts=dfts,log=text)
        output=os.environ.get('OPENEMS_TEST_OUTPUT')
        if output:
            digest=hashlib.sha256(Path(binary).read_bytes()).hexdigest()
            name=engine+('-pec-' if pec else '-surface-')+digest[:12]
            target=Path(output)/name
            target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copytree(folder,target)
            record={k:v for k,v in result.items() if k not in ['log','dfts']}
            record.update(binary_sha256=digest,xml_sha256=hashlib.sha256((target/'model.xml').read_bytes()).hexdigest(),
                engine=engine,pec=pec,dfts={k:[v.real,v.imag] for k,v in dfts.items()})
            (target/'result.json').write_text(json.dumps(record,indent=2)+'\n')
        return result


@unittest.skipUnless(os.environ.get('OPENEMS_BIN'),'Set OPENEMS_BIN for compiled regressions')
class SurfaceStability(unittest.TestCase):
    def test_graded_grid_uses_surface_capacitance_for_timestep(self):
        # Full decay on a grid which diverges with the previous surface update.
        # The default Rennings timestep must include the physical exterior C.
        binary=os.path.abspath(os.environ['OPENEMS_BIN'])
        results=[]
        for engine in ['multithreaded','tiling']:
            result=solve(binary,engine)
            self.assertEqual(result['exit'],0,result['log'][-2000:])
            self.assertIsNone(result['reason'])
            self.assertTrue(result['converged'])
            self.assertFalse(result['capped'])
            self.assertTrue(result['finite'] and result['nonzero'])
            self.assertLess(result['dt'],1.59e-14)
            results.append(result)
        self.assertAlmostEqual(results[0]['dt'],results[1]['dt'],delta=1e-20)
        for name in results[0]['dfts']:
            a,b=[r['dfts'][name] for r in results]
            self.assertLess(abs(a-b)/max(abs(a),abs(b)),3e-5)
        baseline=os.environ.get('OPENEMS_BASELINE_BIN')
        if baseline:
            result=solve(os.path.abspath(baseline))
            self.assertEqual(result['reason'],'energy growth')
            self.assertFalse(result['converged'])

    def test_no_surface_keeps_existing_timestep(self):
        result=solve(os.path.abspath(os.environ['OPENEMS_BIN']),pec=True)
        self.assertEqual(result['exit'],0,result['log'][-2000:])
        self.assertTrue(result['converged'] and result['finite'] and result['nonzero'])
        self.assertIsNone(result['reason'])
        self.assertFalse(result['capped'])
        # Printed log precision; the native unspecialized grid timestep is fixed.
        self.assertAlmostEqual(result['dt'],1.59336e-14,delta=1e-19)


if __name__=='__main__':unittest.main()

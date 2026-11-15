"""A controller-side setup exception must not strand operator workers."""
import subprocess
import sys
import tempfile
import unittest


class OperatorSetupCleanupTest(unittest.TestCase):
    def test_repeated_invalid_surface_setup_returns_instead_of_hanging(self):
        code=r'''
from pathlib import Path
import sys
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
for index in range(2):
    fdtd=openEMS(NrTS=10)
    fdtd.SetGaussExcite(100e9,50e9)
    fdtd.SetBoundaryCond(['PEC']*6)
    csx=ContinuousStructure();fdtd.SetCSX(csx)
    grid=csx.GetGrid();grid.SetDeltaUnit(1e-3)
    for axis in 'xyz':grid.AddLine(axis,np.linspace(0,2,21))
    metal=csx.AddMetal('invalid')
    metal.AddAttribute('ExperimentalOpaqueSurfaceV1','-1;10,1e-12')
    metal.AddBox([.8,.8,.8],[1.2,1.2,1.2])
    try:
        fdtd.Run(str(Path(sys.argv[1])/str(index)),engine='multithreaded',numThreads=2)
    except RuntimeError as error:
        assert 'opaque surface' in str(error).lower(),str(error)
    else:
        raise AssertionError('Invalid passive fit was not rejected')
print('returned from both setup failures')
'''
        with tempfile.TemporaryDirectory() as folder:
            # A broken join used to block indefinitely. Run in a child so a
            # regression fails within a bound instead of hanging the suite.
            process=subprocess.run([sys.executable,'-c',code,folder],
                capture_output=True,text=True,timeout=15)
            self.assertEqual(process.returncode,0,process.stdout+process.stderr)
            self.assertIn('returned from both setup failures',process.stdout)


if __name__=='__main__':unittest.main()

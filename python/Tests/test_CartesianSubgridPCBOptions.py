"""Combined child knots, material sampling and count-budget configuration."""
from pathlib import Path
import os
import tempfile
import unittest
import xml.etree.ElementTree as ET
from CSXCAD import ContinuousStructure
from openEMS import openEMS


class TestPCBOptions(unittest.TestCase):
    def test_joint_native_setup_xml_and_independent_reset(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'combined.xml'
            base = [0, 1, 2, 3, 4]
            children = [0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
            fdtd = openEMS(NrTS=8, EndCriteria=0, CartesianSubgridBox=[1, 1, 1, 3, 3, 3],
                          CartesianSubgridMaxLeaves=4000000, CartesianSubgridMaterialSamples=4,
                          CartesianSubgridChildLines=[children] * 3)
            fdtd.SetGaussExcite(1e9, 1e9)
            fdtd.SetBoundaryCond(['PEC'] * 6)
            csx = ContinuousStructure()
            mesh = csx.GetGrid()
            mesh.SetDeltaUnit(.001)
            for axis in 'xyz':
                mesh.SetLines(axis, base)
            csx.AddMaterial('dielectric', epsilon=3.25, kappa=.02).AddBox([0, 0, 0], [4, 4, 4])
            fdtd.SetCSX(csx)
            self.assertTrue(fdtd.Write2XML(str(path)))
            setup = ET.parse(path).find('FDTD')
            self.assertEqual(setup.get('CartesianSubgridMaxLeaves'), '4000000')
            self.assertEqual(setup.get('CartesianSubgridMaterialSamples'), '4')
            self.assertIsNotNone(setup.find('CartesianSubgridChildLines'))
            restored = openEMS()
            self.assertTrue(restored.ReadFromXML(str(path)))
            cwd = os.getcwd()
            try:
                restored.Run(str(Path(tmp) / 'native'), engine='basic', numThreads=1, exact_endcriteria=True)
            finally:
                os.chdir(cwd)
            restored.SetCartesianSubgridMaxLeaves(1000000)
            self.assertTrue(restored.Write2XML(str(path)))
            setup = ET.parse(path).find('FDTD')
            self.assertIsNone(setup.get('CartesianSubgridMaxLeaves'))
            self.assertEqual(setup.get('CartesianSubgridMaterialSamples'), '4')
            restored.SetCartesianSubgridMaterialSamples(0)
            self.assertTrue(restored.Write2XML(str(path)))
            setup = ET.parse(path).find('FDTD')
            self.assertIsNone(setup.get('CartesianSubgridMaterialSamples'))
            self.assertIsNotNone(setup.find('CartesianSubgridChildLines'))


if __name__ == '__main__':
    unittest.main()

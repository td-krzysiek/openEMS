"""Optional count guard: strict default, atomic rejection and XML reset."""
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS


class TestMaxLeaves(unittest.TestCase):
    def test_roundtrip_reset_and_rejected_update(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'budget.xml'
            fdtd = openEMS(CartesianSubgridBox=[-1, -1, -1, 1, 1, 1], CartesianSubgridMaxLeaves=4000000)
            fdtd.SetGaussExcite(1e9, 1e9)
            fdtd.SetCSX(ContinuousStructure())
            self.assertTrue(fdtd.Write2XML(str(path)))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridMaxLeaves'), '4000000')
            other = openEMS()
            self.assertTrue(other.ReadFromXML(str(path)))
            self.assertTrue(other.Write2XML(str(path)))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridMaxLeaves'), '4000000')
            for value in [True, False, np.bool_(True), 0, -1, 3.5, '4000000', 357913942, 2**64]:
                with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                    other.SetCartesianSubgridMaxLeaves(value)
            self.assertTrue(other.Write2XML(str(path)))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridMaxLeaves'), '4000000')
            other.SetCartesianSubgridMaxLeaves(1000000)
            self.assertTrue(other.Write2XML(str(path)))
            self.assertIsNone(ET.parse(path).find('FDTD').get('CartesianSubgridMaxLeaves'))
            self.assertTrue(fdtd.ReadFromXML(str(path)))
            self.assertTrue(fdtd.Write2XML(str(path)))
            self.assertIsNone(ET.parse(path).find('FDTD').get('CartesianSubgridMaxLeaves'))

    def test_bad_xml(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.xml'
            for value in ['0', '-1', '357913942', '4000000.0', '4000000,2', 'nan', 'inf', '2suffix', '9999999999999999999999999']:
                path.write_text('<openEMS><ContinuousStructure/><FDTD CartesianSubgridMaxLeaves="' + value + '"/></openEMS>')
                with self.subTest(value=value), self.assertRaises((ValueError, RuntimeError)):
                    openEMS().ReadFromXML(str(path))


if __name__ == '__main__':
    unittest.main()

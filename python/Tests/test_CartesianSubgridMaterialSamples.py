"""Explicit sampled-material API/XML, with strict defaults and invalid guards."""
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from openEMS import openEMS
from CSXCAD import ContinuousStructure


class TestMaterialSamples(unittest.TestCase):
    def test_roundtrip_and_clear(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'sampled.xml'
            fdtd = openEMS(CartesianSubgridBox=[-1, -1, -1, 1, 1, 1], CartesianSubgridMaterialSamples=4)
            fdtd.SetGaussExcite(1e9, 1e9)
            fdtd.SetCSX(ContinuousStructure())
            self.assertTrue(fdtd.Write2XML(str(path)))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridMaterialSamples'), '4')
            restored = openEMS()
            self.assertTrue(restored.ReadFromXML(str(path)))
            self.assertTrue(restored.Write2XML(str(path)))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridMaterialSamples'), '4')
            restored.SetCartesianSubgridMaterialSamples(0)
            self.assertTrue(restored.Write2XML(str(path)))
            self.assertIsNone(ET.parse(path).find('FDTD').get('CartesianSubgridMaterialSamples'))

    def test_invalid_python_and_xml(self):
        fdtd = openEMS()
        for value in [True, False, np.bool_(True), -1, 1, 3, 16, float('nan'), float('inf'), '4']:
            with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                fdtd.SetCartesianSubgridMaterialSamples(value)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'invalid.xml'
            for value in ['-1', '1', '3', '4.5', '4suffix', '4,8', 'nan']:
                path.write_text('<openEMS><ContinuousStructure/><FDTD CartesianSubgridMaterialSamples="' + value + '"/></openEMS>')
                with self.subTest(value=value), self.assertRaises((ValueError, RuntimeError)):
                    openEMS().ReadFromXML(str(path))


if __name__ == '__main__':
    unittest.main()

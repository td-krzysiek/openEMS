"""API and XML contract for the experimental planar Cartesian grid."""
import os
import tempfile
import unittest
import xml.etree.ElementTree as ET

from CSXCAD import ContinuousStructure
from openEMS import openEMS


class TestCartesianMultiGrid(unittest.TestCase):
    def test_finite_coordinate_required(self):
        fdtd = openEMS()
        for value in (float('nan'), float('inf'), -float('inf')):
            with self.subTest(value=value), self.assertRaises(ValueError):
                fdtd.SetCartesianMultiGrid(value)

    def test_xml_roundtrip_and_reset(self):
        with tempfile.TemporaryDirectory() as directory:
            selected = os.path.join(directory, 'selected.xml')
            restored = os.path.join(directory, 'restored.xml')
            default = os.path.join(directory, 'default.xml')
            fdtd = openEMS(CartesianMultiGridX=0.0)
            fdtd.SetGaussExcite(1e9, 1e9)
            csx = ContinuousStructure()
            grid = csx.GetGrid()
            for axis in 'xyz':
                grid.SetLines(axis, [-4, -3, -2, -1, 0, 1, 2, 3, 4])
            fdtd.SetCSX(csx)
            self.assertTrue(fdtd.Write2XML(selected))
            self.assertEqual(float(ET.parse(selected).find('FDTD').get('CartesianMultiGridX')), 0.0)
            other = openEMS()
            self.assertTrue(other.ReadFromXML(selected))
            self.assertTrue(other.Write2XML(restored))
            self.assertEqual(float(ET.parse(restored).find('FDTD').get('CartesianMultiGridX')), 0.0)
            plain = openEMS()
            plain.SetGaussExcite(1e9, 1e9)
            plain.SetCSX(ContinuousStructure())
            self.assertTrue(plain.Write2XML(default))
            self.assertTrue(other.ReadFromXML(default))
            self.assertTrue(other.Write2XML(restored))
            self.assertIsNone(ET.parse(restored).find('FDTD').get('CartesianMultiGridX'))

    def test_invalid_xml_is_an_exception(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'bad.xml')
            with open(path, 'w', encoding='utf-8') as stream:
                stream.write('<openEMS><ContinuousStructure/><FDTD CartesianMultiGridX="bad"/></openEMS>')
            with self.assertRaises(ValueError):
                openEMS().ReadFromXML(path)


if __name__ == '__main__':
    unittest.main()

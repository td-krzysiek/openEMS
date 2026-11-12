"""API / optional XML attribute contract for local XYZ refinement."""
import os
import tempfile
import unittest
import xml.etree.ElementTree as ET

from CSXCAD import ContinuousStructure
from openEMS import openEMS


class TestCartesianSubgrid(unittest.TestCase):
    def test_invalid_box(self):
        fdtd = openEMS()
        for box in ([0, 0, 0], [0, 0, 0, 1, 1, float('nan')],
                    [0, 0, 0, 1, 1, float('inf')], [1, 0, 0, 0, 1, 1]):
            with self.subTest(box=box), self.assertRaises(ValueError):
                fdtd.SetCartesianSubgrid(box)

    def test_xml_roundtrip_and_reset(self):
        box = [-4, -4, -4, 4, 4, 4]
        with tempfile.TemporaryDirectory() as directory:
            selected = os.path.join(directory, 'selected.xml')
            restored = os.path.join(directory, 'restored.xml')
            default = os.path.join(directory, 'default.xml')
            fdtd = openEMS(CartesianSubgridBox=box)
            fdtd.SetGaussExcite(1e9, 1e9)
            fdtd.SetCSX(ContinuousStructure())
            self.assertTrue(fdtd.Write2XML(selected))
            other = openEMS()
            self.assertTrue(other.ReadFromXML(selected))
            self.assertTrue(other.Write2XML(restored))
            value = ET.parse(restored).find('FDTD').get('CartesianSubgridBox')
            self.assertEqual([float(v) for v in value.split(',')], box)
            plain = openEMS()
            plain.SetGaussExcite(1e9, 1e9)
            plain.SetCSX(ContinuousStructure())
            self.assertTrue(plain.Write2XML(default))
            self.assertTrue(other.ReadFromXML(default))
            self.assertTrue(other.Write2XML(restored))
            self.assertIsNone(ET.parse(restored).find('FDTD').get('CartesianSubgridBox'))

    def test_invalid_xml_never_drops_bad_tokens(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'bad.xml')
            for box in ('bad,0,0,0,4,4,4', '0,0,0,4,4,4,', '0,0,0,4,4,4bad'):
                with open(path, 'w', encoding='utf-8') as stream:
                    stream.write('<openEMS><ContinuousStructure/><FDTD CartesianSubgridBox="' + box + '"/></openEMS>')
                with self.subTest(box=box), self.assertRaises(ValueError):
                    openEMS().ReadFromXML(path)


if __name__ == '__main__':
    unittest.main()

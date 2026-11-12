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
            fdtd = openEMS(CartesianSubgridBox=box, CartesianSubgridRatio=[1, 2, 8])
            fdtd.SetGaussExcite(1e9, 1e9)
            fdtd.SetCSX(ContinuousStructure())
            self.assertTrue(fdtd.Write2XML(selected))
            other = openEMS()
            self.assertTrue(other.ReadFromXML(selected))
            self.assertTrue(other.Write2XML(restored))
            value = ET.parse(restored).find('FDTD').get('CartesianSubgridBox')
            self.assertEqual([float(v) for v in value.split(',')], box)
            self.assertEqual(ET.parse(restored).find('FDTD').get('CartesianSubgridRatio'), '1,2,8')
            # Legacy files without a ratio must restore the original 2x2x2.
            tree = ET.parse(selected)
            del tree.find('FDTD').attrib['CartesianSubgridRatio']
            tree.write(selected)
            self.assertTrue(other.ReadFromXML(selected))
            self.assertTrue(other.Write2XML(restored))
            self.assertEqual(ET.parse(restored).find('FDTD').get('CartesianSubgridRatio'), '2,2,2')
            plain = openEMS()
            plain.SetGaussExcite(1e9, 1e9)
            plain.SetCSX(ContinuousStructure())
            self.assertTrue(plain.Write2XML(default))
            self.assertTrue(other.ReadFromXML(default))
            self.assertTrue(other.Write2XML(restored))
            self.assertIsNone(ET.parse(restored).find('FDTD').get('CartesianSubgridBox'))
            self.assertIsNone(ET.parse(restored).find('FDTD').get('CartesianSubgridRatio'))

    def test_axis_ratios(self):
        fdtd = openEMS()
        for ratio in ([1, 2], [0, 2, 2], [2, 3, 2], [2, 2.5, 4],
                      [-1, 2, 2], [32, 2, 2], [2, float('nan'), 4],
                      [2, float('inf'), 4], [[2], [2], [2]]):
            with self.subTest(ratio=ratio), self.assertRaises(ValueError):
                fdtd.SetCartesianSubgrid([-4, -4, -4, 4, 4, 4], ratio=ratio)
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'ratio.xml')
            fdtd.SetGaussExcite(1e9, 1e9)
            fdtd.SetCSX(ContinuousStructure())
            fdtd.SetCartesianSubgrid([-4, -4, -4, 4, 4, 4], ratio=[2, 4, 8])
            fdtd.SetCartesianSubgrid([-2, -2, -2, 2, 2, 2])
            self.assertTrue(fdtd.Write2XML(path))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridRatio'), '2,4,8')
            with self.assertRaises(ValueError):
                fdtd.SetCartesianSubgrid([4, -4, -4, -4, 4, 4], ratio=[1, 1, 1])
            self.assertTrue(fdtd.Write2XML(path))
            self.assertEqual(ET.parse(path).find('FDTD').get('CartesianSubgridRatio'), '2,4,8')

    def test_invalid_ratio_xml(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'bad.xml')
            for ratio in ('2,2', '2,2,2,', '2,,2,4', '2.0,2,4', '2,2,4bad',
                          '-1,2,4', '2,3,4', '2,2,9999999999999999999999'):
                with open(path, 'w', encoding='utf-8') as stream:
                    stream.write('<openEMS><ContinuousStructure/><FDTD '
                                 'CartesianSubgridBox="0,0,0,4,4,4" '
                                 'CartesianSubgridRatio="' + ratio + '"/></openEMS>')
                with self.subTest(ratio=ratio), self.assertRaises(ValueError):
                    openEMS().ReadFromXML(path)
            with open(path, 'w', encoding='utf-8') as stream:
                stream.write('<openEMS><ContinuousStructure/><FDTD CartesianSubgridRatio="1,2,4"/></openEMS>')
            with self.assertRaises(ValueError):
                openEMS().ReadFromXML(path)

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

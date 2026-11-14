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

    def test_multiple_regions_roundtrip_and_updates(self):
        boxes = [[-6, -4, -4, -2, 4, 4], [2, -4, -4, 6, 4, 4]]
        ratios = [[2, 2, 4], [4, 1, 2]]
        fdtd = openEMS(CartesianSubgrids=[{'box': boxes[0], 'ratio': ratios[0]}])
        fdtd.AddCartesianSubgrid(boxes[1], ratio=ratios[1])
        fdtd.SetGaussExcite(1e9, 1e9)
        fdtd.SetCSX(ContinuousStructure())
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'multi.xml')
            def saved_regions():
                self.assertTrue(fdtd.Write2XML(path))
                opts = ET.parse(path).find('FDTD')
                self.assertIsNone(opts.get('CartesianSubgridBox'))
                self.assertIsNone(opts.get('CartesianSubgridRatio'))
                return [dict(child.attrib) for child in opts.findall('CartesianSubgrids/Region')]
            original = saved_regions()
            self.assertEqual(len(original), 2)
            for i, region in enumerate(original):
                self.assertEqual([float(v) for v in region['Box'].split(',')], boxes[i])
                self.assertEqual([int(v) for v in region['Ratio'].split(',')], ratios[i])
            other = openEMS()
            self.assertTrue(other.ReadFromXML(path))
            restored = os.path.join(directory, 'restored.xml')
            self.assertTrue(other.Write2XML(restored))
            self.assertEqual([dict(e.attrib) for e in ET.parse(restored).findall('FDTD/CartesianSubgrids/Region')], original)
            for box in (boxes[0], [-4, -2, -2, -3, 2, 2], [-2, -4, -4, 0, 4, 4]):
                with self.subTest(box=box), self.assertRaises(ValueError):
                    fdtd.AddCartesianSubgrid(box)
                self.assertEqual(saved_regions(), original)
            with self.assertRaises(ValueError):
                fdtd.SetCartesianSubgridRatio([2, 2, 2])
            self.assertEqual(saved_regions(), original)
            with self.assertRaises(ValueError):
                fdtd.SetCartesianSubgrids([{'box': boxes[0]}, {'box': boxes[0]}])
            self.assertEqual(saved_regions(), original)
            # The singular setter explicitly replaces the entire region list.
            fdtd.SetCartesianSubgrid([-4, -4, -4, 4, 4, 4], ratio=[1, 2, 8])
            self.assertTrue(fdtd.Write2XML(path))
            opts = ET.parse(path).find('FDTD')
            self.assertIsNone(opts.find('CartesianSubgrids'))
            self.assertEqual(opts.get('CartesianSubgridRatio'), '1,2,8')
            self.assertTrue(other.ReadFromXML(path))
            self.assertTrue(other.Write2XML(restored))
            self.assertIsNone(ET.parse(restored).find('FDTD/CartesianSubgrids'))

    def test_invalid_multiple_regions(self):
        fdtd = openEMS()
        for regions in ([], [{}], [{'box': [0,0,0,2,2,2], 'ratios': [2,2,2]}],
                        [{'box': [0,0,0,2,2,2], 'ratio': [2,3,2]}],
                        [{'box': [0,0,0,2,2,2]}] * 65):
            with self.subTest(regions=regions), self.assertRaises(ValueError):
                fdtd.SetCartesianSubgrids(regions)
        with self.assertRaises(ValueError):
            openEMS(CartesianSubgrids=[{'box': [0,0,0,2,2,2]}], CartesianSubgridBox=[0,0,0,2,2,2])

    def test_invalid_multi_xml(self):
        valid = '<Region Box="1,1,1,3,3,3" Ratio="2,2,4"/>'
        invalid = ('<CartesianSubgrids/>',
                   '<CartesianSubgrids><Unknown/></CartesianSubgrids>',
                   '<CartesianSubgrids><Region Ratio="2,2,2"/></CartesianSubgrids>',
                   '<CartesianSubgrids><Region Box="1,1,1,3,3,3" Ratoi="2,4,8"/></CartesianSubgrids>',
                   '<CartesianSubgrids><Region Box="1,1,1,3,3,3bad"/></CartesianSubgrids>',
                   '<CartesianSubgrids><Region Box="1,1,1,3,3,3" Ratio="2,,4"/></CartesianSubgrids>',
                   '<CartesianSubgrids>' + valid + valid + '</CartesianSubgrids>',
                   '<CartesianSubgrids>' + valid + '</CartesianSubgrids><CartesianSubgrids>' + valid + '</CartesianSubgrids>')
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'bad.xml')
            for contents in invalid:
                with open(path, 'w', encoding='utf-8') as stream:
                    stream.write('<openEMS><ContinuousStructure/><FDTD>' + contents + '</FDTD></openEMS>')
                with self.subTest(contents=contents), self.assertRaises(ValueError):
                    openEMS().ReadFromXML(path)
            with open(path, 'w', encoding='utf-8') as stream:
                stream.write('<openEMS><ContinuousStructure/><FDTD CartesianSubgridBox="1,1,1,3,3,3">'
                             '<CartesianSubgrids>' + valid + '</CartesianSubgrids></FDTD></openEMS>')
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

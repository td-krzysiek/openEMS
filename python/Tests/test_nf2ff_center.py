"""NF2FF phase-centre units and closed-box preflight, without a simulation."""
import unittest
from types import SimpleNamespace
import numpy as np
from openEMS.nf2ff import nf2ff


def recording():
    obj=object.__new__(nf2ff)
    grid=SimpleNamespace(GetMeshType=lambda:0,GetDeltaUnit=lambda:1e-3)
    obj.CSX=SimpleNamespace(GetGrid=lambda:grid)
    obj.start=[-5,-5,-3];obj.stop=[5,5,3]
    obj.directions=[True]*6;obj.mirror=[0]*6
    return obj


class NF2FFCenterTest(unittest.TestCase):
    def test_metres_independent_of_mm_drawing_unit(self):
        obj=recording()
        np.testing.assert_allclose(obj._validate_center([0,0,.00013589]),[0,0,.00013589])
        with self.assertRaisesRegex(ValueError,'metres'):
            obj.CalcNF2FF('unused',24e9,[0],[0],center=[0,0,.13589])

    def test_invalid_vector_and_boundary(self):
        obj=recording()
        for center in [[0,0],[0,0,np.nan],[0,0,.003]]:
            with self.assertRaises(ValueError):obj._validate_center(center)

    def test_open_mirrored_and_cylindrical_boxes_keep_their_conventions(self):
        for kind in ['open','mirror','cylinder']:
            obj=recording()
            if kind=='open':obj.directions[0]=False
            elif kind=='mirror':obj.mirror[0]=1
            else:obj.CSX.GetGrid().GetMeshType=lambda:1
            np.testing.assert_allclose(obj._validate_center([0,0,.003]),[0,0,.003])


if __name__=='__main__':unittest.main()

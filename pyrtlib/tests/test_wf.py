from unittest import TestCase

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal
from pyrtlib.weighting_functions import WeightingFunctions
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.utils import ppmv2gkg, mr2rh


class Test(TestCase):
    def wf_computation(self, frq: np.ndarray):
        z, p, d, t, md = atmp.gl_atm(atmp.TROPICAL)
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100
        
        wf = WeightingFunctions(z, p, t, rh)
        wf.frequencies = frq

        return wf.generate_wf(), z

    def test_f_89(self):
        wgt, z = self.wf_computation(np.array([89.0]))
        np.testing.assert_almost_equal(z[np.argmax(wgt)], 0., decimal=15)
        
    def test_f_57(self):
        wgt, z = self.wf_computation(np.array([57.290344]))
        np.testing.assert_almost_equal(z[np.argmax(wgt)], 18., decimal=15)
        
    def test_f_57_6(self):
        wgt, z = self.wf_computation(np.array([57.660544]))
        np.testing.assert_almost_equal(z[np.argmax(wgt)], 27.5, decimal=15)

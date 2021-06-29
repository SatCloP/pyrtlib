import os
from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_equal
from pyrtlib.absmod_uncertainty import absmod_uncertainties_perturb
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):
    def test_absmod_uncertainties_perturb_min(self):
        amu = absmod_uncertainties_perturb(['gamma_a'], 'min', 1)
        expected0 = 0.016501363999999998
        assert_equal(expected0, amu['gamma_a'].uncer[0])
        expected1 = 0.015001240000000
        assert_equal(expected1, amu['gamma_a'].uncer[1])

        expected0 = 2.698723076000000
        assert_equal(expected0, amu['gamma_a'].value[0])
        expected1 = 2.929742172000000
        assert_equal(expected1, amu['gamma_a'].value[1])

    def test_absmod_uncertainties_perturb_max(self):
        amu = absmod_uncertainties_perturb(['gamma_a'], 'max', 0)
        expected0 = 2.7152244399999996
        assert_equal(expected0, amu['gamma_a'].value[0])

        amu = absmod_uncertainties_perturb(['gamma_a'], 'max', 1)
        expected1 = 2.9597446520000004
        assert_equal(expected1, amu['gamma_a'].value[1])

    def test_pyrtlib_uncertainty_gamma_a(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        frq = np.arange(20, 201, 1)

        amu = absmod_uncertainties_perturb(['gamma_a'], 'min', index=1)

        rte = BTCloudRTE(z, p, t, rh, frq, amu=amu)
        rte.init_absmdl('uncertainty')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tbtotal_uncertainty_gamma_a_min.csv"))
        assert_allclose(df.tbtotal.values, df_expected.gamma_a)

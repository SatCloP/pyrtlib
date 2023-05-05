import os
from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_equal
from pyrtlib.uncertainty import AbsModUncertainty, SpectroscopicParameter
from pyrtlib.uncertainty import covariance_matrix
from pyrtlib.atmospheric_profiles import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):
    def test_absmod_uncertainties_perturb_min(self):
        amu = AbsModUncertainty.parameters_perturbation(['gamma_a'], 'min', 1)
        expected0 = 0.016501363999999998
        assert_equal(expected0, amu['gamma_a'].uncer[0])
        expected1 = 0.015001240000000
        assert_equal(expected1, amu['gamma_a'].uncer[1])

        expected0 = 2.698723076000000
        assert_equal(expected0, amu['gamma_a'].value[0])
        expected1 = 2.929742172000000
        assert_equal(expected1, amu['gamma_a'].value[1])

    def test_absmod_uncertainties_perturb_max(self):
        amu = AbsModUncertainty.parameters_perturbation(['gamma_a'], 'max', 0)
        expected0 = 2.7152244399999996
        assert_equal(expected0, amu['gamma_a'].value[0])

        amu = AbsModUncertainty.parameters_perturbation(['gamma_a'], 'max', 1)
        expected1 = 2.9597446520000004
        assert_equal(expected1, amu['gamma_a'].value[1])

    def test_pyrtlib_uncertainty_gamma_a(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        frq = np.arange(20, 201, 1)

        amu = AbsModUncertainty.parameters_perturbation(
            ['gamma_a'], 'min', index=1)

        rte = TbCloudRTE(z, p, t, rh, frq, amu=amu)
        rte.init_absmdl('uncertainty')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tbtotal_uncertainty_gamma_a_min.csv"))
        assert_allclose(df.tbtotal.values, df_expected.gamma_a)

    def test_cov_p(self):
        n = 106
        r111 = covariance_matrix.R17_111
        r112 = covariance_matrix.R17_112

        assert_allclose(r111[0:n, 0:n], r112[0:n, 0:n])

    def test_spectroscopic_params(self):
        parameters = SpectroscopicParameter.parameters()
        v = parameters['w2a'].value
        assert_equal(v, 1.2)

    def test_add_spectr_params(self):
        parameters = SpectroscopicParameter.parameters()
        parameters['test'] = SpectroscopicParameter(
            2.3, 0.001, 'unitless', 'Tretyakov, JMS, 2016')
        parameters['test'].value
        assert_equal(parameters['test'].value, 2.3)

    def test_edit_spectroscopic_params(self):
        parameters = SpectroscopicParameter.parameters()
        parameters['w2a'].value = 1.4
        assert_equal(parameters['w2a'].value, 1.4)

    def test_uncertainty_propagation(self):
        parameters = SpectroscopicParameter.parameters()
        u = AbsModUncertainty.uncertainty_propagation(
            "A/B", parameters['delta_a'].value[0],
            parameters['gamma_a'].value[0],
            parameters['delta_a'].uncer[0],
            parameters['gamma_a'].uncer[0])
        expected = (-0.012229016120066702, 0.001854236371206569)
        assert_equal(expected, u)

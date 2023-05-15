import os
from unittest import TestCase
import pytest
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_equal
from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, O3AbsModel
from pyrtlib.uncertainty import AbsModUncertainty, SpectroscopicParameter
from pyrtlib.uncertainty import covariance_matrix
from pyrtlib.atmospheric_profiles import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh, ppmv_to_moleculesm3

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class Test(TestCase):
    def test_absmod_uncertainties_perturb_min(self):
        SpectroscopicParameter._initialize()
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
        SpectroscopicParameter._initialize()
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

        water_sp = SpectroscopicParameter.water_parameters("R17")
        oxygen_sp = SpectroscopicParameter.oxygen_parameters("R18")

        parameters = {**water_sp, **oxygen_sp}
        SpectroscopicParameter.set_parameters(parameters)

        amu = AbsModUncertainty.parameters_perturbation(
            ['gamma_a'], 'min', index=1)

        rte = TbCloudRTE(z, p, t, rh, frq, amu=amu)
        rte.init_absmdl('R17')
        O2AbsModel.model = 'R18'
        O2AbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tbtotal_uncertainty_gamma_a_min.csv"))
        assert_allclose(df.tbtotal.values, df_expected.gamma_a)

    def test_pyrtlib_uncertainty_wO3(self):
        z, p, _, t, md = atmp.gl_atm(atmp.US_STANDARD)

        o3n_ppmv = md[:, atmp.O3]
        o3n = np.zeros(z.shape)
        for k in range(0, len(z)):
            o3n[k] = ppmv_to_moleculesm3(o3n_ppmv[k], p[k] * 100.0, t[k])

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        frq = np.arange(20, 201, 1)

        water_sp = SpectroscopicParameter.water_parameters("R21SD")
        oxygen_sp = SpectroscopicParameter.oxygen_parameters("R20")
        ozono_sp = SpectroscopicParameter.ozono_parameters("R18")

        parameters = {**water_sp, **oxygen_sp, **ozono_sp}
        SpectroscopicParameter.set_parameters(parameters)

        amu = AbsModUncertainty.parameters_perturbation(
            ['O3_X'], 'min', index=0)

        rte = TbCloudRTE(z, p, t, rh, frq, o3n=o3n, amu=amu)
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        O3AbsModel.model = 'R18'
        O3AbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tbtotal_uncertainty_gamma_a_min.csv"))
        assert_allclose(df.tbtotal.values, df_expected.o3_x)

    def test_cov_p(self):
        n = 106
        r111 = covariance_matrix.R17_111
        r112 = covariance_matrix.R17_112

        assert_allclose(r111[0:n, 0:n], r112[0:n, 0:n])

    # @pytest.mark.skip(reason="skipping")
    def test_spectroscopic_params(self):
        oxygen_sp = SpectroscopicParameter.oxygen_parameters("R18")
        v = oxygen_sp['w2a'].value
        assert_equal(v, 1.2)

    def test_add_spectr_params(self):
        parameters = SpectroscopicParameter.water_parameters("R17")
        parameters['test'] = SpectroscopicParameter(
            2.3, 0.001, 'unitless', 'Tretyakov, JMS, 2016')
        parameters['test'].value
        assert_equal(parameters['test'].value, 2.3)

    # @pytest.mark.skip(reason="skipping")
    def test_edit_spectroscopic_params(self):
        parameters = SpectroscopicParameter.oxygen_parameters("R18")
        parameters['w2a'].value = 1.4
        assert_equal(parameters['w2a'].value, 1.4)

    # @pytest.mark.skip(reason="skipping")
    def test_set_parameters(self):
        parameters = SpectroscopicParameter.water_parameters("R17")
        parameters['gamma_a'].value[0] = 2.688
        parameters['gamma_a'].uncer[0] = 0.039
        # SpectroscopicParameter.set_parameters(parameters)
        assert_equal(parameters['gamma_a'].value[0], 2.688)

    def test_uncertainty_propagation(self):
        parameters = SpectroscopicParameter.water_parameters("R19")
        parameters['delta_a'].uncer = np.array([0.005])
        parameters['gamma_a'].uncer = np.array([0.022])
        u = AbsModUncertainty.uncertainty_propagation(
            "A/B", parameters['delta_a'].value[0],
            parameters['gamma_a'].value[0],
            parameters['delta_a'].uncer[0],
            parameters['gamma_a'].uncer[0])
        expected = (-0.0122267506483882929, 0.0018552168412088639)
        assert_equal(expected, u)
        u = AbsModUncertainty.uncertainty_propagation(
            "aA", parameters['delta_a'].value[0],
            parameters['gamma_a'].value[0],
            parameters['delta_a'].uncer[0],
            parameters['gamma_a'].uncer[0])
        expected = (-0.033, 0.005)
        assert_equal(expected, u)
        u = AbsModUncertainty.uncertainty_propagation(
            "aA+bB", parameters['delta_a'].value[0],
            parameters['gamma_a'].value[0],
            parameters['delta_a'].uncer[0],
            parameters['gamma_a'].uncer[0])
        expected = (2.666, 0.02256102834535695)
        assert_equal(expected, u)


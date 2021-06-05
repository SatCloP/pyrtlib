import os
from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
from pyrtlib.absmod_uncertainty import absmod_uncertainties_perturb, AMU
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):
    def test_absmod_uncertainties_perturb_min(self):
        amu = absmod_uncertainties_perturb(['gamma_a'], 'min', 2)
        expected0 = 0.016501364000000
        assert_allclose(expected0, amu['gamma_a'].uncer[0])
        expected1 = 0.015001240000000
        assert_allclose(expected1, amu['gamma_a'].uncer[1])

        expected0 = 2.698723076000000
        assert_allclose(expected0, amu['gamma_a'].value[0])
        expected1 = 2.929742172000000
        assert_allclose(expected1, amu['gamma_a'].value[1])

    def test_absmod_uncertainties_perturb_max(self):
        amu = absmod_uncertainties_perturb(['gamma_a'], 'max', 1)
        expected0 = 2.715224440000000
        assert_allclose(expected0, amu['gamma_a'].value[0])

        amu = absmod_uncertainties_perturb(['gamma_a'], 'max', 2)
        expected1 = 2.959744652000000
        assert_allclose(expected1, amu['gamma_a'].value[1])

    def test_pyrtlib_uncertainty_gamma_a(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        amu = absmod_uncertainties_perturb(['gamma_a'], 'min', index=1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
        rte.init_absmdl('uncertainty')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tbtotal_uncertainty_gamma_a_min.csv"))
        assert_allclose(df.tbtotal, df_expected.gamma_a, atol=0)

    def foo():
        # O2 parameters
        O2_parameters = {'O2S': [],
                         'X05': [],
                         'WB300': [],
                         'O2gamma': [],
                         'Y300': [],
                         'O2_V': []}

        for i in range(0, 34):
            O2_parameters['O2gamma'].append(i)
            O2_parameters['Y300'].append(i)
            O2_parameters['O2_V'].append(i)

        uncertainties = {'O2S': [np.mean(AMU['O2S'].uncer / AMU['O2S'].value) * 100],
                         'X05': [AMU['X05'].uncer],
                         'WB300': [AMU['WB300'].uncer],
                         'O2gamma': AMU['O2gamma'].uncer[0: 34].tolist(),
                         'Y300': AMU['Y300'].uncer[0: 34].tolist(),
                         'O2_V': AMU['O2_V'].uncer[0: 34].tolist()}
        uncertainties_all = [item for sublist in list(uncertainties.values()) for item in sublist]

        # H2O parameters
        HO2_parameters = {'con_Cf_factr': [],
                         'con_Cs_factr': [],
                         'gamma_a': [0],
                         'S': [0],
                         'con_Xf': [],
                         'SR': [0],
                         'con_Xs': []}

        uncertainties = {'con_Cf_factr': [AMU['con_Cf'].value * AMU['con_Cf_factr'].uncer],
                         'con_Cs_factr': [AMU['con_Cs'].value * AMU['con_Cs_factr'].uncer],
                         'gamma_a': [AMU['gamma_a'].uncer[0]],
                         'S': [AMU['S'].uncer[0]],
                         'con_Xf': [AMU['con_Xf'].uncer],
                         'SR': [AMU['SR'].uncer[0]],
                         'con_Xs': [AMU['con_Xs'].uncer]}
        uncertainties_all = [item for sublist in list(uncertainties.values()) for item in sublist]

        for k, v in O2_parameters.items():
            if v:
                for i in v:
                    amu_p = absmod_uncertainties_perturb([k], 'min', index=i)
            else:
                amu_p = absmod_uncertainties_perturb([k], 'min')


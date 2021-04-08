import os
# from pathlib import Path
from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import tb_cloud_rte
from pyrtlib.utils import ppmv2gkg, mr2rh

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_tb_cloud_with_rose19sd(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        denliq = np.zeros(z.shape)
        denice = np.zeros(z.shape)
        cldh = np.zeros((2, 0))

        df = tb_cloud_rte(z, p, t, rh, denliq, denice, cldh, frq, ang,
                          absmdl='rose19sd',
                          ray_tracing=True,
                          from_sat=False)
        df = df.set_index(frq)
        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ground_ros03_19sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    # @pytest.mark.skip(reason="rose03 not completly implemented yet")
    def test_tb_cloud_with_rose03(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        denliq = np.zeros(z.shape)
        denice = np.zeros(z.shape)
        cldh = np.zeros((2, 0))

        df = tb_cloud_rte(z, p, t, rh, denliq, denice, cldh, frq, ang,
                          absmdl='rose03',
                          ray_tracing=True,
                          from_sat=False)
        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ground_ros03_19sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros03, atol=0)

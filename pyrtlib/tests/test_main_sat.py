import os
from pyrtlib.absmodel import H2OAbsModel, O2AbsModel
# from pathlib import Path
from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh, import_lineshape

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose19sd_atm(self):
        d = {'tropical': atmp.TROPICAL,
             'midlat_summer': atmp.MIDLATITUDE_SUMMER,
             'midlat_winter': atmp.MIDLATITUDE_WINTER,
             'subarctic_summer': atmp.SUBARCTIC_SUMMER,
             'subarctic_winter': atmp.SUBARCTIC_WINTER,
             'us_standard': atmp.US_STANDARD}

        for k, v in d.items():
            z, p, _, t, md = atmp.gl_atm(v)

            gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
            rh = mr2rh(p, t, gkg)[0] / 100

            ang = np.array([90.])
            frq = np.arange(20, 201, 1)

            rte = BTCloudRTE(z, p, t, rh, frq, ang)
            rte.init_absmdl('rose19sd')
            df = rte.execute()

            df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tbtotal_atm_rose19sd.csv"))
            assert_allclose(df.tbtotal, df_expected[k], atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose19sd(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose19sd')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose19(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose19')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose16(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose16')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros16, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_makarov11(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose16')
        O2AbsModel.model = 'makarov11'
        O2AbsModel.o2ll = import_lineshape('o2ll_{}'.format(O2AbsModel.model))
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros16, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    # @pytest.mark.skip(reason="rose03 not completly implemented yet")
    def test_pyrtlib_sat_rose03(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose03')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros03, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose17(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose17')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros17, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose20(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose20')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros20, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose20sd(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose20sd')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros20sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose18(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose18')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros18, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose98(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose98')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros98, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose98_cloudy(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        denliq = np.zeros(z.shape)
        denice = np.zeros(z.shape)
        cldh = np.zeros((2, 2))
        # build a cloud
        ib = 1
        it = 3
        denliq[ib:it + 1] = 10 * np.ones((it - ib + 1))
        cldh[:, 0] = np.array([z[ib], z[it]])
        ib = 29
        it = 31
        denice[ib:it + 1] = 0.1 * np.ones((it - ib + 1))
        cldh[:, 1] = np.array([z[ib], z[it]])

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('rose98')
        rte.cloudy = True
        rte.init_cloudy(cldh, denice, denliq)
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "cloudy_tb_tot_ros98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros98, atol=0)

    def test_pyrtlib_sat_no_raytracing(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.ray_tracing = False
        rte.init_absmdl('rose19sd')

        df = rte.execute()
        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19sd, atol=0)
import os
# from pathlib import Path
from unittest import TestCase
import pytest
import numpy as np
import pandas as pd
from datetime import datetime
from numpy.testing import assert_allclose
from pyrtlib.absmodel import H2OAbsModel, LiqAbsModel, O2AbsModel, O3AbsModel
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.apiwebservices import WyomingUpperAir, ERA5Reanalysis, IGRAUpperAir
from pyrtlib.utils import ppmv2gkg, ppmv_to_moleculesm3, mr2rh, import_lineshape, dewpoint2rh, kgkg_to_kgm3, pressure_to_height, to_kelvin

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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.mak11, atol=0)

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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros20sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_rose21sd(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq)
        rte.init_absmdl('rose20')
        H2OAbsModel.model = 'rose21sd'
        H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros21sd, atol=0)

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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
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

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.rosen, atol=0)

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
        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19sd, atol=0)

    def test_pyrtlib_sat_rose21_wO3(self):
        z, p, _, t, md = atmp.gl_atm(atmp.US_STANDARD)

        o3n_ppmv = md[:, atmp.O3]
        o3n = np.zeros(z.shape)
        for k in range(0, len(z)):
            o3n[k] = ppmv_to_moleculesm3(o3n_ppmv[k], p[k] * 100.0, t[k])
        # o3n = np.interp(z,o3n,z)
        # o3n_matlab = np.loadtxt("/Users/slarosa/Downloads/o3.csv", delimiter=',')
        # assert_allclose(o3n, o3n_matlab, atol=0)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang, o3n)
        rte.init_absmdl('rose20')
        H2OAbsModel.model = 'rose21sd'
        H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
        O3AbsModel.model = 'rose18'
        O3AbsModel.o3ll = import_lineshape('o3ll_{}'.format(O3AbsModel.model))
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros21sd_wo3.csv"))
        assert_allclose(df.tbtotal, df_expected.ros21sd_wo3, atol=0)

    def test_pyrtlib_sat_rose21sd_wyoming_es(self):
        date = datetime(2021, 4, 22, 12)
        station = 'LIRE'
        df_w = WyomingUpperAir.request_data(date, station)

        z, p, t, gkg = df_w.height.values / 1000, \
                    df_w.pressure.values, \
                    to_kelvin(df_w.temperature.values), \
                    df_w.mixr.values

        rh = dewpoint2rh(df_w.dewpoint, df_w.temperature).values

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.emissivity = 0.6
        rte.init_absmdl('rose20')
        H2OAbsModel.model = 'rose21sd'
        H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal_wyoming, atol=0)
        
    def test_pyrtlib_sat_rose21sd_igra2_es(self):
        date = datetime(2020, 6, 1, 12)
        station = 'SPM00008221'
        df_igra2, header = IGRAUpperAir.request_data(date, station)
        
        df_igra2 = df_igra2[df_igra2.pressure.notna() & 
                            df_igra2.temperature.notna() & 
                            df_igra2.dewpoint.notna() & 
                            df_igra2.height.notna()]

        z, p, t = df_igra2.height.values / 1000, df_igra2.pressure.values, to_kelvin(df_igra2.temperature.values)
        
        rh = dewpoint2rh(df_igra2.dewpoint, df_igra2.temperature).values

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.emissivity = 0.6
        rte.init_absmdl('rose20')
        H2OAbsModel.model = 'rose21sd'
        H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal_igra2, atol=0)
        
    def test_pyrtlib_sat_rose21sd_igra2_beg2021_es(self):
        date = datetime(2021, 10, 2, 0)
        station = 'ASM00094610'
        df_igra2, header = IGRAUpperAir.request_data(date, station, beg2021=True)
        
        df_igra2 = df_igra2[df_igra2.pressure.notna() & 
                            df_igra2.temperature.notna() & 
                            df_igra2.dewpoint.notna() & 
                            df_igra2.height.notna()]

        z, p, t = df_igra2.height.values / 1000, df_igra2.pressure.values, to_kelvin(df_igra2.temperature.values)
        
        rh = dewpoint2rh(df_igra2.dewpoint, df_igra2.temperature).values

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = BTCloudRTE(z, p, t, rh, frq, ang)
        rte.emissivity = 0.6
        rte.init_absmdl('rose20')
        H2OAbsModel.model = 'rose21sd'
        H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal_igra2_beg2021, atol=0)

    def test_pyrtlib_sat_rose21sd_ERA5_cloudy(self):
        lonlat = (15.8158, 38.2663)
        nc_file = os.path.join(THIS_DIR, "data", "era5_reanalysis-2019-06-25T12:00:00.nc")
        df_era5 = ERA5Reanalysis.read_data(nc_file, lonlat)

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        cldh = np.empty((2, 1))
        cldh[:, 0] = np.array([np.min(df_era5.z), np.max(df_era5.z)])

        total_mass = 1 - df_era5.ciwc.values - df_era5.clwc.values - df_era5.crwc.values - df_era5.cswc.values
        denice = df_era5.ciwc.values * (1/total_mass) * kgkg_to_kgm3(df_era5.q.values * (1/total_mass),
                                                    df_era5.p.values, df_era5.t.values) * 1000
        denliq = df_era5.clwc.values * (1/total_mass) * kgkg_to_kgm3(df_era5.q.values * (1/total_mass),
                                                    df_era5.p.values, df_era5.t.values) * 1000


        rte = BTCloudRTE(df_era5.z.values, df_era5.p.values, df_era5.t.values, df_era5.rh.values, frq, ang)
        rte.init_absmdl('rose20')
        rte.init_cloudy(cldh, denice, denliq)
        H2OAbsModel.model = 'rose21sd'
        H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
    
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_ERA5_cloudy.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal, atol=0)
        
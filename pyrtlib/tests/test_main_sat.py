import os
# from pathlib import Path
from unittest import TestCase
import pytest
import numpy as np
import pandas as pd
from datetime import datetime
from numpy.testing import assert_allclose
from pyrtlib.absorption_model import H2OAbsModel, LiqAbsModel, O2AbsModel, O3AbsModel
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.apiwebservices import WyomingUpperAir, ERA5Reanalysis, IGRAUpperAir
from pyrtlib.utils import ppmv2gkg, ppmv_to_moleculesm3, mr2rh, import_lineshape, dewpoint2rh, kgkg_to_kgm3, pressure_to_height, to_kelvin

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_R19SD_atm(self):
        d = atmp.atm_profiles()

        for k, v in d.items():
            z, p, _, t, md = atmp.gl_atm(k)

            gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
            rh = mr2rh(p, t, gkg)[0] / 100

            ang = np.array([90.])
            frq = np.arange(20, 201, 1)

            rte = TbCloudRTE(z, p, t, rh, frq, ang)
            rte.init_absmdl('R19SD')
            df = rte.execute()

            df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tbtotal_atm_rose19sd.csv"))
            assert_allclose(df.tbtotal, df_expected[v.lower().replace(" ", "_")], atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R19SD(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R19SD')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R19(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R19')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R16(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R16')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros16, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_makarov11(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R16')
        O2AbsModel.model = 'makarov11'
        O2AbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.mak11, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R03(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R03')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros03, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R17(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R17')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros17, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R20(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R20')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros20, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R20SD(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R20SD')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros20sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    @pytest.mark.skip(reason="skipping")
    def test_pyrtlib_sat_R21SD(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq)
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros21sd, atol=0)

    def test_pyrtlib_sat_R18(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R18')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.ros18, atol=0)

    def test_pyrtlib_sat_R98(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R98')
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros03_16_17_18_19_19sd_20_20sd_98_mak11_21sd.csv"))
        assert_allclose(df.tbtotal, df_expected.rosen, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_sat_R98_cloudy(self):
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

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.init_absmdl('R98')
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

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.ray_tracing = False
        rte.init_absmdl('R19SD')

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

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang, o3n)
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        O3AbsModel.model = 'R18'
        O3AbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_ros21sd_wo3.csv"))
        assert_allclose(df.tbtotal, df_expected.ros21sd_wo3, atol=0)

    def test_pyrtlib_sat_R21SD_wyoming_es(self):
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

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.emissivity = 0.6
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal_wyoming, atol=0)
        
    def test_pyrtlib_sat_R21SD_igra2_es(self):
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

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.emissivity = 0.6
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal_igra2, atol=0)
        
    def test_pyrtlib_sat_R21SD_igra2_beg2021_es(self):
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

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.emissivity = 0.6
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal_igra2_beg2021, atol=0)

    def test_pyrtlib_sat_R21SD_ERA5_cloudy(self):
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


        rte = TbCloudRTE(df_era5.z.values, df_era5.p.values, df_era5.t.values, df_era5.rh.values, frq, ang)
        rte.init_absmdl('R20')
        rte.init_cloudy(cldh, denice, denliq)
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
    
        df = rte.execute()

        df_expected = pd.read_csv(
            os.path.join(THIS_DIR, "data", "tb_tot_rose21sd_RAOB_es.csv"))
        assert_allclose(df.tbtotal, df_expected.tbtotal, atol=0)
        
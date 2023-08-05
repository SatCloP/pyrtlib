import os
# from pathlib import Path
from unittest import TestCase

import numpy as np
import pandas as pd
from numpy.testing import assert_allclose, assert_equal
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.absorption_model import H2OAbsModel
from pyrtlib.apiwebservices import ERA5Reanalysis
from pyrtlib.utils import ppmv2gkg, mr2rh, import_lineshape

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_pyrtlib_ground_R19SD(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq, ang)
        rte.satellite = False
        rte.init_absmdl('R19SD')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ground_ros03_19sd_21sd_era5.csv"))
        assert_allclose(df.tbtotal, df_expected.ros19sd, atol=0)

    # @pytest.mark.datafiles(DATA_DIR)
    # @pytest.mark.skip(reason="R03 not completly implemented yet")
    def test_pyrtlib_ground_R03(self):
        z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(z, p, t, rh, frq)
        rte.satellite = False
        rte.init_absmdl('R03')
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ground_ros03_19sd_21sd_era5.csv"))
        assert_allclose(df.tbtotal.values, df_expected.ros03)

    def test_pyrtlib_ground_R21SD_ERA5(self):
        lonlat = (15.8158, 38.2663)
        nc_file = os.path.join(THIS_DIR, "data", "era5_reanalysis-2019-06-25T12:00:00.nc")
        df_era5 = ERA5Reanalysis.read_data(nc_file, lonlat)

        ang = np.array([90.])
        frq = np.arange(20, 201, 1)

        rte = TbCloudRTE(df_era5.z.values, df_era5.p.values, df_era5.t.values, df_era5.rh.values, frq, ang)
        rte.satellite = False
        rte.init_absmdl('R20')
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        df = rte.execute()

        df_expected = pd.read_csv(os.path.join(THIS_DIR, "data", "tb_tot_ground_ros03_19sd_21sd_era5.csv"))
        assert_allclose(df.tbtotal, df_expected.rose21sd_era5, atol=0)

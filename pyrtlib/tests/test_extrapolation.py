import os
# from pathlib import Path
from unittest import TestCase
import pytest
import numpy as np
import pandas as pd
from datetime import datetime
from numpy.testing import assert_allclose
from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, O3AbsModel
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.climatology import ProfileExtrapolation
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.apiwebservices import WyomingUpperAir, ERA5Reanalysis, IGRAUpperAir
from pyrtlib.utils import ppmv2gkg, ppmv_to_moleculesm3, mr2rh, dewpoint2rh, kgkg_to_kgm3, to_kelvin

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class Test(TestCase):
    def test_profile_extr(self):
        date = datetime(2020, 6, 1, 12)
        station = 'SPM00008221'
        df_igra2, header = IGRAUpperAir.request_data(date, station)

        df_igra2 = df_igra2[df_igra2.pressure.notna() &
                            df_igra2.temperature.notna() &
                            df_igra2.dewpoint.notna() &
                            df_igra2.height.notna()]

        z, p, t = df_igra2.height.values / 1000, df_igra2.pressure.values, to_kelvin(df_igra2.temperature.values)
        rh = dewpoint2rh(df_igra2.dewpoint, df_igra2.temperature).values
        assert len(z) < 25
        assert min(p) > 10
        ex = ProfileExtrapolation()
        zz, pp, tt, rhh = ex.profile_extrapolation(header.latitude.values[0], 6, z, (p, t, rh))
        assert len(zz) > 25
        assert min(pp) < 10
import os
# from pathlib import Path
from unittest import TestCase
from datetime import datetime
from numpy.testing import assert_allclose
from pyrtlib.climatology import ProfileExtrapolation
from pyrtlib.apiwebservices import IGRAUpperAir
from pyrtlib.utils import dewpoint2rh, to_kelvin

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

        z, p, t = df_igra2.height.values / \
            1000, df_igra2.pressure.values, to_kelvin(
                df_igra2.temperature.values)
        rh = dewpoint2rh(df_igra2.dewpoint, df_igra2.temperature).values
        assert len(z) < 25
        assert min(p) > 10
        ex = ProfileExtrapolation()
        zz, pp, _, _ = ex.profile_extrapolation(
            header.latitude.values[0], 6, z, (p, t, rh))
        assert len(zz) > 25
        assert min(pp) < 10
        
    def test_standard_temp(self):
        ex = ProfileExtrapolation()
        t = ex.standard_temperature(12, 280)
        assert_allclose(t, 216.65)
        
    def test_standard_pres(self):
        ex = ProfileExtrapolation()
        p = ex.standard_pressure(12, 280, 1013)
        assert_allclose(p, 193.99616049)

    def test_standard_wv(self):
        ex = ProfileExtrapolation()
        wv = ex.standard_water_vapour_pressure(12)
        assert_allclose(wv, 0.018586351836920856)
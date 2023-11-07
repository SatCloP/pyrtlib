import os
import numpy as np
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
        
    def test_standard_wvd(self):
        ex = ProfileExtrapolation()
        wvd = ex.standard_water_vapour_density(12, 1, 7.5)
        assert_allclose(wvd, 4.608159264996157e-05)
        
    def test_pressure(self):
        ex = ProfileExtrapolation()
        p = ex.pressure(38, np.array([1, 2, 5]), 'winter')
        assert_allclose(p, np.array([899.398 , 789.5947, 518.1532]))
        
        ex = ProfileExtrapolation()
        p = ex.pressure(38, np.array([1, 2, 5]), 'summer')
        assert_allclose(p, np.array([905.1263, 805.1632, 551.6491]))
        
    def test_temperature(self):
        ex = ProfileExtrapolation()
        t = ex.temperature(38, np.array([1, 2, 5]), 'winter')
        assert_allclose(t, np.array([268.9265, 264.7771, 250.2181]))
        
        t = ex.temperature(38, np.array([1, 2, 5]), 'summer')
        assert_allclose(t, np.array([289.69681, 284.26764, 267.12705]))
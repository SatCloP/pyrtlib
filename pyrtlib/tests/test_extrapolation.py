import os
# from pathlib import Path
from unittest import TestCase
from datetime import datetime
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
        zz, pp, tt, rhh = ex.profile_extrapolation(
            header.latitude.values[0], 6, z, (p, t, rh))
        assert len(zz) > 25
        assert min(pp) < 10

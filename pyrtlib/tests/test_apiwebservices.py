import os
# from pathlib import Path
from unittest import TestCase
import pytest
import numpy as np
import pandas as pd
from datetime import datetime
from numpy.testing import assert_allclose, assert_equal
from pyrtlib.apiwebservices import WyomingUpperAir, ERA5Reanalysis, IGRAUpperAir

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_wyom_get_stations(self):
        df = WyomingUpperAir.get_stations()
        station_id = np.array(['03005', '03238', '03808', '03882', '03918', '06011', '06458',
                               '10035', '10113', '10184', '10393', '10410', '10548', '10739',
                               '10868', '11520', '11747', '11952', '12120', '12374', '12425',
                               '12843', '12982', '13275', '13388', '14240', '14430', '15420',
                               '16045', '16064', '16144', '16245', '16332', '16429', '16546',
                               '16622', '16754', '17030', '17064', '17196', '17220', '17240',
                               '17351', '22820', '22845', '26038', '26075', '26298', '26477',
                               '26629', '26781', '27038', '27199', '27459', '27594', '27713',
                               '27707', '27730', '27962', '27995', '34009', '34247', '34467',
                               '37011', '60390', '60715'])

        assert_equal(df.station_id.values, station_id)

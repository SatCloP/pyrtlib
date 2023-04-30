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
        df = WyomingUpperAir.get_stations('pac')
        b = df[df.station_id=='45004'].station_name == 'Kings Park'

        assert b.values

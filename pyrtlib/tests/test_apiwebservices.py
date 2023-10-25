import os
# from pathlib import Path
from unittest import TestCase
from pyrtlib.apiwebservices import WyomingUpperAir

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_wyom_get_stations(self):
        df = WyomingUpperAir.get_stations('pac')
        b = df[df.station_id=='45004'].station_name == 'Kings Park'

        assert b.values

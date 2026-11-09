import os
# from pathlib import Path
from unittest import TestCase
from datetime import datetime
import tempfile
import pandas as pd
import pytest
from pyrtlib.apiwebservices import WyomingUpperAir
from pyrtlib.apiwebservices.erafive import ERA5Reanalysis

# TEST_DIR = Path(__file__).parent
# DATA_DIR = os.path.join(TEST_DIR, 'data')
THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):

    # @pytest.mark.datafiles(DATA_DIR)
    def test_wyom_get_stations(self):
        df = WyomingUpperAir.get_stations()
        assert 'station_id' in df.columns
        assert 'station_name' in df.columns
        
        df_station = WyomingUpperAir.find_station('Africa')
        print(df_station)
        assert not df_station.empty
        
    @pytest.mark.skip()
    def test_get_era5(self):
        lonlat = (15.13, 37.87)
        date = datetime(2022, 7, 22, 12)
        nc_file = ERA5Reanalysis.request_data(tempfile.gettempdir(), date, lonlat)
        print(nc_file)
        df_era5 = ERA5Reanalysis.read_data(nc_file, lonlat)
        
        assert df_era5.attrs['units']['t'] == 'K'
        assert df_era5.attrs['units']['p'] == 'hPa'
        assert df_era5.attrs['units']['rh'] == '%'

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

    def test_wyom_find_station_not_found(self):
        time = datetime(2026, 11, 7, 12)

        df = WyomingUpperAir.find_station(
            'StationThatDoesNotExist',
            time=time
        )

        assert isinstance(df, pd.DataFrame)
        assert df.empty

    def test_wyom_find_station_case_insensitive(self):
        time = datetime(2026, 11, 7, 12)

        df_upper = WyomingUpperAir.find_station('AFRICA', time=time)
        df_lower = WyomingUpperAir.find_station('africa', time=time)

        assert not df_upper.empty
        assert not df_lower.empty

        assert df_upper['station_id'].tolist() == \
            df_lower['station_id'].tolist()

    def test_wyom_get_stations_without_datetime(self):
        df = WyomingUpperAir.get_stations()

        assert isinstance(df, pd.DataFrame)
        assert 'station_id' in df.columns
        assert 'station_name' in df.columns
        assert not df.empty
        
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

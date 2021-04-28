"""
Read and Download data from ERA5Reanalysis Reanalysis model data.
To download ERA5 dataset it is necessary to configure a API key.
Step-by-step ti create a API key can be found to https://cds.climate.copernicus.eu/api-how-to
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

import os
import math
from datetime import datetime
from typing import Optional, Tuple

import cdsapi
import numpy as np
import pandas as pd
from netCDF4 import Dataset

from ..utils import pressure_to_height


class ERA5Reanalysis:
    """Read and Download data from ERA5 CDS Reanalysis model data"""

    @classmethod
    def read_data(cls, file: str, lonlat: tuple) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, datetime]:
        """Read data from the ERA5Reanalysis dataset.

        Args:
            file (str): The netcdf file
            lonlat (tuple): longitude and latitude

        Returns:
            [type]: [description]
        """
        ERAFIVE = cls()
        nc = Dataset(file)
        lats = nc.variables['latitude'][:]
        lons = nc.variables['longitude'][:]
        idx_lat = ERAFIVE.find_nearest(lats, lonlat[1])
        idx_lon = ERAFIVE.find_nearest(lons, lonlat[0])
        pres = np.asarray(nc.variables['level'][:])
        temp = np.asarray(nc.variables['t'][:, :, idx_lat, idx_lon])
        rh = np.asarray(nc.variables['r'][:, :, idx_lat, idx_lon])
        z = pressure_to_height(pres) / 1000
        date = pd.to_datetime(nc.variables['time'][:], origin='1900-01-01 00:00:00.0', unit='h')

        return np.flip(z), np.flip(pres), np.flip(temp[0]), np.flip(rh[0]), date

    def find_nearest(self, array, value):
        """Find index of nearest coordinate

        Args:
            array ([type]): [description]
            value ([type]): [description]

        Returns:
            [type]: [description]
        """
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            return idx - 1
        else:
            return idx

    @staticmethod
    def request_data(path: str, time: datetime, lonlat: tuple, offset: Optional[np.float] = 0.3) -> str:
        """Download ERA5Reanalysis data from the Copernicus Climate Change Service.

        Args:
            path (str): The output directory
            time (datetime): The date and time of the desired observation.
            lonlat (tuple): The coordinatre in degrees, longitude and latitude
            offset (Optional[np.float], optional): The offset to apply to coordinates to get the extent. Defaults to 0.3.

        Returns:
            str: The path to downloaded netcdf file
        """
        # North, West, South, Est
        extent = [lonlat[1] + offset, lonlat[0] - offset, lonlat[1] - offset, lonlat[0] + offset]
        nc_file_name = 'era5_reanlysis-{}.nc'.format(time.isoformat())
        nc_file = os.path.join(path, nc_file_name)
        c = cdsapi.Client()
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'variable': [
                    'relative_humidity', 'temperature',
                ],
                'pressure_level': [
                    '1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150',
                    '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600',
                    '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950',
                    '975', '1000',
                ],
                'year': time.year,
                'month': time.month,
                'day': time.day,
                'time': '{}:00'.format(time.hour),
                'area': extent,
                'grid': [0.25, 0.25],
                'format': 'netcdf',
            },
            nc_file)

        return nc_file

"""
Read and Download data from ERA5Reanalysis Reanalysis model data.
To download ERA5 dataset it is necessary to configure an API key.
Step-by-step to create an API key can be found to https://cds.climate.copernicus.eu/api-how-to
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

import os
import warnings
from datetime import datetime
from typing import Optional

try:
    import cdsapi
except ModuleNotFoundError as e:
    warnings.warn(
        "Module CDSAPI must be installed to download ERA5 Reanalysis dataset.")
import numpy as np
# from scipy.spatial import cKDTree
from sklearn.neighbors import BallTree
import pandas as pd
from netCDF4 import Dataset

from ..utils import atmospheric_tickness


class ERA5Reanalysis:
    """Read and Download data from ERA5 CDS Reanalysis model data"""

    @classmethod
    def read_data(cls, file: str, lonlat: tuple) -> pd.DataFrame:
        """Read data from the ERA5 Reanalysis dataset.

        Args:
            file (str): The netcdf file
            lonlat (tuple): longitude and latitude

        Returns:
            pandas.DataFrame: Dataframe containing the variables retrieved.

        .. note:: Variables name and units information are reported within the attribute `units` of
            the returned dataframe (see example below).
        
        Example:
            .. code-block:: python

                >>> from pyrtlib.apiwebservices import ERA5Reanalysis
                >>> lonlat = (15.8158, 38.2663)
                >>> date = datetime(2020, 2, 22, 12)
                >>> nc_file = ERA5Reanalysis.request_data(tempfile.gettempdir(), date, lonlat)
                >>> df_era5 = ERA5Reanalysis.read_data(nc_file, lonlat)
                >>> df_era5.attrs['units']
                {'p': 'hPa',
                'z': 'km',
                't': 'K',
                'rh': '%',
                'clwc': 'kg kg-1',
                'ciwc': 'kg kg-1',
                'crwc': 'kg kg-1',
                'cswc': 'kg kg-1',
                'o3': 'kg kg-1',
                'q': 'kg kg-1'}

        .. note:: To convert specific cloud water content (CLWC) or specific cloud ice water content (CIWC) 
            from kg kg-1 to g m-3 using this function :py:meth:`pyrtlib.utils.kgkg_to_gm3`
        """
        ERAFIVE = cls()
        nc = Dataset(file)
        lats = nc.variables['latitude'][:]
        lons = nc.variables['longitude'][:]
        # idx_lat = ERAFIVE.find_nearest(lats, lonlat[1])
        # idx_lon = ERAFIVE.find_nearest(lons, lonlat[0])
        idx, _ = ERAFIVE._find_nearest(lons, lats, lonlat)

        pres = np.asarray(nc.variables['level'][:])
        temp = np.asarray(nc.variables['t'][:, :, idx, idx])
        # RH in decimal
        rh = np.asarray(nc.variables['r'][:, :, idx, idx]) / 100
        clwc = np.asarray(nc.variables['clwc'][:, :, idx, idx])
        ciwc = np.asarray(nc.variables['ciwc'][:, :, idx, idx])
        crwc = np.asarray(nc.variables['crwc'][:, :, idx, idx])
        cswc = np.asarray(nc.variables['cswc'][:, :, idx, idx])
        ozone = np.asarray(nc.variables['o3'][:, :, idx, idx])
        q = np.asarray(nc.variables['q'][:, :, idx, idx])

        z = atmospheric_tickness(np.flip(pres), np.flip(
            temp[0]), np.flip(q[0]))  # Altitude in km

        date = pd.to_datetime(
            nc.variables['time'][:], origin='1900-01-01 00:00:00.0', unit='h')

        df = pd.DataFrame({'p': np.flip(pres),
                           'z': z,
                           't': np.flip(temp[0]),
                           'rh': np.flip(rh[0]),
                           'clwc': np.flip(clwc[0]),
                           'ciwc': np.flip(ciwc[0]),
                           'crwc': np.flip(crwc[0]),
                           'cswc': np.flip(cswc[0]),
                           'o3': np.flip(ozone[0]),
                           'q': np.flip(q[0]),
                           'time': np.repeat(date, len(pres))
                           })
        df.attrs['units'] = {'p': 'hPa',
                             'z': 'km',
                             't': 'K',
                             'rh': '%',
                             'clwc': 'kg kg-1',
                             'ciwc': 'kg kg-1',
                             'crwc': 'kg kg-1',
                             'cswc': 'kg kg-1',
                             'o3': 'kg kg-1',
                             'q': 'kg kg-1'}

        return df  # , (idx, dist), np.stack((lons, lats))

    def _find_nearest(self, lons, lats, point):
        """Find index of nearest coordinate

        Args:
            array ([type]): [description]
            value ([type]): [description]

        Returns:
            [type]: [description]
        """
        # idx = np.searchsorted(array, value, side="left")
        # if np.amax(array) < value < np.amin(array):
        #     warnings.warn('Out of boundingbox: value {} is out of dataset extent'.format(value))
        # if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
        #     return idx - 1
        # else:
        #     return idx

        # scipy
        # combined_x_y_arrays = np.dstack([lats.ravel(), lons.ravel()])[0]
        # # points = list(point.transpose())
        #
        # mytree = cKDTree(combined_x_y_arrays)
        # dist, indexes = mytree.query(point)

        ball = BallTree(np.radians(
            np.stack((lats, lons), axis=-1)), metric='haversine')
        dist, indexes = ball.query(np.radians(
            np.array(np.flip(point)).reshape(1, 2)), k=1)

        return indexes[0][0], dist

    @staticmethod
    def request_data(path: str, time: datetime, lonlat: tuple, resolution: Optional[float] = 0.25, offset: Optional[float] = 0.4) -> str:
        """Download ERA5Reanalysis data from the Copernicus Climate Change Service.

        Args:
            path (str): The output directory
            time (datetime): The date and time of the desired observation.
            lonlat (tuple): The coordinatre in degrees, longitude and latitude
            resolution (Optional[float], optional): The pixel size of the requested grid data. Defaults to 0.25.
            offset (Optional[float], optional): The offset to apply to coordinates to get the extent. Defaults to 0.3.

        Returns:
            str: The path to downloaded netcdf file
        """
        # North, West, South, Est
        extent = [lonlat[1] + offset, lonlat[0] - offset,
                  lonlat[1] - offset, lonlat[0] + offset]
        nc_file_name = 'era5_reanalysis-{}.nc'.format(time.isoformat())
        nc_file = os.path.join(path, nc_file_name)

        variables = ['relative_humidity', 'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content',
                     'specific_humidity', 'specific_rain_water_content', 'specific_snow_water_content', 'ozone_mass_mixing_ratio', 'temperature']
        c = cdsapi.Client()
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'variable': variables,
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
                'grid': [resolution, resolution],
                'format': 'netcdf',
            },
            nc_file)

        return nc_file

"""
Read and Download data from EPS Datastore by EUMETSAT.
To download products from EUMETSAT it is necessary to configure a API key.
To get cinsumer_key and consumer_secret you have to registered to EPS system.
Credentials may be found to the following link https://api.eumetsat.int/api-key/
"""

__author__ = ''
__date__ = 'March 2023'
__copyright__ = '(C) 2021, CNR-IMAA'

import os
from pathlib import Path
import warnings
from datetime import datetime
from typing import Optional, Tuple, Any

try:
    import eumdac
except ModuleNotFoundError as e:
    warnings.warn("Module EUMDAC must be installed to download EPS dataset.")
import numpy as np
# from scipy.spatial import cKDTree
from sklearn.neighbors import BallTree
import pandas as pd
from netCDF4 import Dataset

from ..utils import pressure_to_height


class EUMETSATProduct:
    """Read and Download product from EUMETSAT Data Store Web Services"""

    def __init__(self) -> None:
        super().__init__()
        self.token = None
        self.datastore = None
        self.datatailoe = None
        self.tailor_models = None

    def authentication(self) -> None:
        rc_file = os.path.join(Path.home(), '.eumdacrc')
        if not os.path.exists(rc_file):
            return warnings.warn("eumdacrc file not found in home directory.")
        with open(rc_file, 'r') as f:
            lines = f.readlines()
            consumer_key = lines[0].split(":")[1].strip()
            consumer_secret = lines[1].split(":")[1].strip()

        credentials = (consumer_key, consumer_secret)
        self.token = eumdac.AccessToken(credentials)
        self.datastore = eumdac.DataStore(self.token)
        self.datatailoe = eumdac.DataTailor(self.token)
        self.tailor_models = eumdac.tailor_models
        print(f"This token '{self.token}' expires {self.token.expiration}")

    def list_collections(self, keyword: str = None) -> Tuple:
        """_summary_

        Args:
            keyword (str, optional): _description_. Defaults to None.

        Returns:
            Tuple: _description_
        """
        collections = self.datastore.collections
        ids = []
        if keyword:
            collections_filtered = []
            for cl in collections:
                if keyword in cl.title:
                    ids.append(cl.metadata['properties']['identifier'])
                    collections_filtered.append(cl)
            collections = collections_filtered
        else:
            for cl in collections:
                ids.append(cl.metadata['properties']['identifier'])

        return (ids, collections)

    def product_search(self, collection: eumdac.collection.Collection, **query: Any) -> eumdac.collection.SearchResults:
        """_summary_

        Args:
            collection (eumdac.collection.Collection): _description_

        Returns:
            eumdac.collection.SearchResults: _description_
        """
        products = collection.search(**query)

        return products

    def custom_chain(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        pass

    @classmethod
    def read_data(cls, file: str, lonlat: tuple) -> pd.DataFrame:
        """Read data from the ERA5 Reanalysis dataset.

        Args:
            file (str): The netcdf file
            lonlat (tuple): longitude and latitude

        Returns:
            pandas.DataFrame: [description]

        .. note:: To convert specific cloud water content (CLWC) or specific cloud ice water content (CIWC) 
            from kg kg-1 to g m-3 using this function :py:meth:`pyrtlib.utils.kgkg_to_gm3`
        """
        ERAFIVE = cls()
        nc = Dataset(file)
        lats = nc.variables['latitude'][:]
        lons = nc.variables['longitude'][:]
        # idx_lat = ERAFIVE.find_nearest(lats, lonlat[1])
        # idx_lon = ERAFIVE.find_nearest(lons, lonlat[0])
        idx, dist = ERAFIVE.find_nearest(lons, lats, lonlat)

        pres = np.asarray(nc.variables['level'][:])
        temp = np.asarray(nc.variables['t'][:, :, idx, idx])
        # RH in decimal
        rh = np.asarray(nc.variables['r'][:, :, idx, idx]) / 100
        clwc = np.asarray(nc.variables['clwc'][:, :, idx, idx])
        ciwc = np.asarray(nc.variables['ciwc'][:, :, idx, idx])
        crwc = np.asarray(nc.variables['crwc'][:, :, idx, idx])
        cswc = np.asarray(nc.variables['cswc'][:, :, idx, idx])
        q = np.asarray(nc.variables['q'][:, :, idx, idx])

        z = pressure_to_height(pres) / 1000  # Altitude in km
        date = pd.to_datetime(
            nc.variables['time'][:], origin='1900-01-01 00:00:00.0', unit='h')

        df = pd.DataFrame({'p': np.flip(pres),
                           'z': np.flip(z),
                           't': np.flip(temp[0]),
                           'rh': np.flip(rh[0]),
                           'clwc': np.flip(clwc[0]),
                           'ciwc': np.flip(ciwc[0]),
                           'crwc': np.flip(crwc[0]),
                           'cswc': np.flip(cswc[0]),
                           'q': np.flip(q[0]),
                           'time': np.repeat(date, len(z))
                           })

        return df  # , (idx, dist), np.stack((lons, lats))

    def find_nearest(self, lons, lats, point):
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
        # extent = [lonlat[1] + offset, lonlat[0] - offset,
        #           lonlat[1] - offset, lonlat[0] + offset]
        nc_file_name = 'era5_reanalysis-{}.nc'.format(time.isoformat())
        nc_file = os.path.join(path, nc_file_name)

        # variables = ['relative_humidity', 'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content',
        #              'specific_humidity', 'specific_rain_water_content', 'specific_snow_water_content', 'temperature']
        return nc_file

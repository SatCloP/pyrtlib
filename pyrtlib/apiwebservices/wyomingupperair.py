"""Read upper air data from the Wyoming archives.
The orignal code of this class can be found in https://github.com/Unidata/siphon
Terms and conditions are as for siphon library license (https://github.com/Unidata/siphon/blob/master/LICENSE)
"""

import warnings
from datetime import datetime
from io import StringIO
from typing import Union, Optional
import requests
from urllib.parse import quote

import pandas as pd
from ..apiwebservices.webservices import HTTPEndPoint

warnings.filterwarnings('ignore', 'Pandas doesn\'t allow columns to be created', UserWarning)


class WyomingUpperAir(HTTPEndPoint):
    """Download and parse data from the University of Wyoming's upper air archive."""

    def __init__(self) -> None:
        """Set up endpoint."""
        super().__init__('https://weather.arcc.uwyo.edu/wsgi/sounding')

    @classmethod
    def request_data(cls, time: datetime, site_id: Union[str, int], **kwargs) -> pd.DataFrame:
        """Retrieve upper air observations from the Wyoming archive.

        Args:
            time (datetime.datetime): The date and time of the desired observation.
            site_id (Union[str, int]): The station identifier for which data should be
                downloaded.

        Returns:
            pandas.DataFrame:  A dataframe containing the data

        .. note:: Variables name and units information are reported within the attribute `units` of
            the returned dataframe (see example below).

        Example:
            .. code-block:: python

                >>> from pyrtlib.apiwebservices import WyomingUpperAir
                >>> from datetime import datetime
                >>> date = datetime(2022, 6, 22, 12)
                >>> station = 'LIRE'
                >>> df = WyomingUpperAir.request_data(date, station)
                >>> df.attrs['units']
                {'pressure': 'hPa',
                 'height': 'meter',
                 'temperature': 'degC',
                 'dewpoint': 'degC',
                 'rh': '%',
                 'mixr': 'g/kg',
                 'station': None,
                 'station_number': None,
                 'time': None,
                 'latitude': 'degrees',
                 'longitude': 'degrees',
                 'elevation': 'meter'}
        """

        endpoint = cls()
        src = kwargs.get('src', 'FM35')
        df = endpoint._get_data(time, site_id, src=src)
        return df

    @classmethod
    def get_stations(cls, time: Optional[datetime] = None) -> pd.DataFrame:
        """Retrieve available stations from the Wyoming archive.

        Args:
            time (Optional[datetime]): Requested sounding date and time.
                If provided, retrieve stations available for that time.

        Returns:
            pandas.DataFrame: A DataFrame of station IDs and names.
        """

        url = 'https://weather.uwyo.edu/cgi-bin/sounding_json'

        params = {}
        if time is not None:
            params['datetime'] = time.strftime('%Y-%m-%d %H:%M:%S')

        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()

        data = resp.json()
        stations = data['stations']

        df = pd.DataFrame(stations)

        df['name'] = df['name'].replace('', 'UNKNOWN').fillna('UNKNOWN')

        df = df.rename(columns={
            'stationid': 'station_id',
            'name': 'station_name'
        })

        return df

    @classmethod
    def find_station(
        cls,
        name: str,
        time: Optional[datetime] = None
    ) -> pd.DataFrame:
        """Search for stations by name (case-insensitive).

        Args:
            name (str): The station name (or part of it) to search for.
            time (Optional[datetime]): Requested sounding date and time.
                If provided, search within stations available for that time.

        Returns:
            pandas.DataFrame: A DataFrame of matching stations.
        """ 

        df = cls.get_stations(time=time)

        mask = df['station_name'].str.contains(
            name, case=False, na=False
        )

        return df[mask].reset_index(drop=True)

    def _get_data(self, time: datetime, site_id: Union[str, int], src: str = 'FM35') -> pd.DataFrame:
        r"""Download and parse upper air observations from an online archive.

        Parameters
        ----------
        time : datetime
            The date and time of the desired observation.

        site_id : str
            The station identifier for which data should be downloaded.

        src : str
            The data source format (e.g. ``FM35`` or ``BUFR``).

        Returns
        -------
            df :  pandas.DataFrame
                containing the data
        """

        raw = self._get_data_raw(time, site_id, src=src)
        df = pd.read_csv(StringIO(raw))

        rename_map = {
            'pressure_hPa': 'pressure',
            'geopotential height_m': 'height',
            'temperature_C': 'temperature',
            'dew point temperature_C': 'dewpoint',
            'ice point temperature_C': 'ice_point_temperature',
            'relative humidity_%': 'rh',
            'humidity wrt ice_%': 'humidity_wrt_ice',
            'mixing ratio_g/kg': 'mixr',
            'wind direction_degree': 'direction',
            'wind speed_m/s': 'speed',
        }
        df = df.rename(columns=rename_map)

        if 'pressure' not in df.columns and 'pressure_hPa' in df.columns:
            df = df.rename(columns={'pressure_hPa': 'pressure'})

        missing = [c for c in ('pressure', 'height', 'temperature', 'mixr') if c not in df.columns]
        if missing:
            raise KeyError(
                f"Missing columns: {missing}. "
                f"Available columns: {list(df.columns)}"
            )

        df['station'] = site_id
        df['time'] = time
        df['title'] = f"{site_id} Observations at {time.strftime('%d %b %Y, %HZ')}"

        df.attrs['units'] = {'pressure': 'hPa',
                             'height': 'meter',
                             'temperature': 'degC',
                             'dewpoint': 'degC',
                             'rh': '%',
                             'mixr': 'g/kg',
                             'direction': 'degrees',
                             'speed': 'm/s',
                             'station': None,
                             'time': None,
                             'latitude': 'degrees',
                             'longitude': 'degrees',
                             'title': None}
        
        return df

    def _get_data_raw(self, time: datetime, site_id: Union[str, int], src: str = 'FM35') -> str:
        """Download data from the University of Wyoming's upper air archive.

        Parameters
        ----------
        time : datetime
            Date and time for which data should be downloaded
        site_id : str
            Site id for which data should be downloaded
        src : str
            The data source format (e.g. ``FM35`` or ``BUFR``).

        Returns
        -------
        text of the server response

        """

        formatted_time = time.strftime('%Y-%m-%d %H:%M:%S')
        path = (f'?datetime={quote(formatted_time)}'
                f'&id={site_id}'
                f'&type=TEXT:CSV'
                f'&src={src}')
        resp = self._get_path(path)
        if '<HTML>' in resp.text[:200] or 'Other data for Station' in resp.text:
            raise ValueError(
                'No data available for {time:%Y-%m-%d %HZ} '
                'for station {stid} (src={src}).'.format(time=time, stid=site_id, src=src))
        return resp.text
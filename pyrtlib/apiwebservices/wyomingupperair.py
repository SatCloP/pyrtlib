"""Read upper air data from the Wyoming archives.
The orignal code of this class can be found in https://github.com/Unidata/siphon
Terms and conditions are as for siphon library license (https://github.com/Unidata/siphon/blob/master/LICENSE)
"""

import warnings
from datetime import datetime
from io import StringIO
from typing import Union, Optional

import pandas as pd
try:
    from bs4 import BeautifulSoup
except ModuleNotFoundError as e:
    warnings.warn("Module BS4 must be installed to download Wyoming RAOB dataset.")
from ..apiwebservices.webservices import HTTPEndPoint

warnings.filterwarnings('ignore', 'Pandas doesn\'t allow columns to be created', UserWarning)


class WyomingUpperAir(HTTPEndPoint):
    """Download and parse data from the University of Wyoming's upper air archive."""

    def __init__(self) -> None:
        """Set up endpoint."""
        super().__init__('http://weather.uwyo.edu/cgi-bin/sounding')

    @classmethod
    def request_data(cls, time: datetime, site_id: Union[str, int], **kwargs) -> pd.DataFrame:
        """Retrieve upper air observations from the Wyoming archive. 

        Args:
            time (datetime.datetime): The date and time of the desired observation.
            site_id (Union[str, int]): The three letter ICAO identifier of the station for which data should be
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
        df = endpoint._get_data(time, site_id)
        return df

    @classmethod
    def get_stations(cls, region: Optional[str] = 'europe') -> pd.DataFrame:
        """Retrieve list of available stations from the Wyoming archive.

        +---------+----------------+
        | region  | name           |
        +---------+----------------+
        | naconf  | North America  |
        +---------+----------------+
        | samer   | South America  |
        +---------+----------------+
        | pac     | South Pacific  |
        +---------+----------------+
        | nz      | New Zealand    |
        +---------+----------------+
        | ant     | Antarctica     |
        +---------+----------------+
        | np      | Arctic         |
        +---------+----------------+
        | europe  | Europe         |
        +---------+----------------+
        | africa  | Africa         |
        +---------+----------------+
        | seasia  | Southeast Asia |
        +---------+----------------+
        | mideast | Mideast        |
        +---------+----------------+

        Args:
            region (Optional[str], optional): The name of region from which to get stations list. Defaults to 'europe'.

        Returns:
            pandas.DataFrame: A dDataFrame of stations id and name
        """

        endpoint = cls()
        endpoint._base = 'http://weather.uwyo.edu/upperair/'

        resp = endpoint._get_path('{}.html'.format(region))
        soup = BeautifulSoup(resp.text, 'html.parser')
        area_tag = soup.find_all('area')
        stations = {}
        for t in area_tag:
            title = t.get('title')
            if title:
                cod, name = title.split(' ')[0], ' '.join(title.split(' ')[1:]).strip()
                stations[cod] = name

        df = pd.DataFrame(list(stations.items()), columns=['station_id', 'station_name'])
        return df

    def _get_data(self, time: datetime, site_id: Union[str, int]) -> pd.DataFrame:
        r"""Download and parse upper air observations from an online archive.

        Parameters
        ----------
        time : datetime
            The date and time of the desired observation.

        site_id : str
            The three letter ICAO identifier of the station for which data should be
            downloaded.

        Returns
        -------
            df :  pandas.DataFrame
                containing the data
        """

        raw_data = self._get_data_raw(time, site_id)
        soup = BeautifulSoup(raw_data, 'html.parser')
        tabular_data = StringIO(soup.find_all('pre')[0].contents[0])
        title = soup.find_all('h2')[0].contents[0]
        col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'rh', 'mixr']
        # df = pd.read_fwf(tabular_data, skiprows=5, usecols=[0, 1, 2, 3, 4, 5], names=col_names, infer_rows=300)
        # change to read_csv function as some data column (e.g. height) couls be corrupted
        # To be further investigated in order to catch the issue
        df = pd.read_csv(tabular_data, header=None, skiprows=5, delim_whitespace=True, usecols=[0, 1, 2, 3, 4, 5], names=col_names)
        # Drop any rows with all NaN values for T, Td, winds
        df = df.dropna(subset=('temperature', 'dewpoint'), how='any').reset_index(drop=True)

        # Parse metadata
        meta_data = soup.find_all('pre')[1].contents[0]
        lines = meta_data.splitlines()

        # If the station doesn't have a name identified we need to insert a
        # record showing this for parsing to proceed.
        if 'Station number' in lines[1]:
            lines.insert(1, 'Station identifier: ')

        station = lines[1].split(':')[1].strip()
        station_number = int(lines[2].split(':')[1].strip())
        sounding_time = datetime.strptime(lines[3].split(':')[1].strip(), '%y%m%d/%H%M')
        latitude = float(lines[4].split(':')[1].strip())
        longitude = float(lines[5].split(':')[1].strip())
        elevation = float(lines[6].split(':')[1].strip())

        df['station'] = station
        df['station_number'] = station_number
        df['time'] = sounding_time
        df['latitude'] = latitude
        df['longitude'] = longitude
        df['elevation'] = elevation
        df['title'] = title

        # Add unit dictionary
        df.attrs['units'] = {'pressure': 'hPa',
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
        return df

    def _get_data_raw(self, time: datetime, site_id: Union[str, int]) -> str:
        """Download data from the University of Wyoming's upper air archive.

        Parameters
        ----------
        time : datetime
            Date and time for which data should be downloaded
        site_id : str
            Site id for which data should be downloaded

        Returns
        -------
        text of the server response

        """
        path = ('?region=naconf&TYPE=TEXT%3ALIST'
                '&YEAR={time:%Y}&MONTH={time:%m}&FROM={time:%d%H}&TO={time:%d%H}'
                '&STNM={stid}').format(time=time, stid=site_id)

        resp = self._get_path(path)
        # See if the return is valid, but has no data
        if resp.text.find('Can\'t') != -1:
            raise ValueError(
                'No data available for {time:%Y-%m-%d %HZ} '
                'for station {stid}.'.format(time=time, stid=site_id))

        return resp.text

# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from typing import Tuple, Optional
import numpy as np

from pyrtlib.utils import rho2rh


class ProfileExtrapolation:

    # This is an abstract class that contains an instance to a modality of the
    # extrapolation.

    def __init__(self, mode: Optional[str] = 'ITU-Annex1'):
        """Initialize what extrapolation to be used.
        Only ITU-835-6 Annex1 has been implemented rigth now.

        Args:
            mode (str, optional): Recommendation of extrapolation. Defaults to 'ITU-Annex1'.

        Raises:
            ValueError: Raise an error if mode is not set or if unsupported.

        References
        ----------
            .. [1] https://www.itu.int/rec/R-REC-P.835/en
        """
        if mode == 'ITU-Annex1':
            self.instance = _ITU835_6()
            self._height = np.hstack((np.linspace(0, 25, 26),
                                      np.linspace(27.5, 50, 10),
                                      np.linspace(55, 100, 10)))
        else:
            raise ValueError("Modality of extrapolation is mandatory")

    @property
    def height(self) -> np.ndarray:
        """Getter/Setter for altitude vector"""
        return self._height

    @height.setter
    def height(self, height: np.ndarray) -> None:
        """Setter for absorption model"""
        if height is not None and isinstance(height, np.ndarray):
            self._height = height
        else:
            raise ValueError("Please enter a valid array for height")

    def temperature(self, lat: float, h: np.ndarray, season: Optional[str] = 'summer') -> np.ndarray:
        """ Determine the temperature at a given latitude and height.

        Method to determine the temperature as a function of altitude and latitude,
        for calculating gaseous attenuation along an Earth-space path. This method
        is recommended when more reliable local data are not available.


        Parameters
        ----------
        lat : numpy.ndarray
            Latitude (degree)
        h : numpy.ndarray
            Height (km)
        season : str
            Season of the year (available values, 'summer', and 'winter').
            Default 'summer'


        Returns
        -------
        T: numpy.ndarray
            Temperature (K)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en

        """
        return self.instance.temperature(lat, h, season)

    def pressure(self, lat: float, h: np.ndarray, season: Optional[str] = 'summer') -> np.ndarray:
        """ Determine the atmospheric pressure at a given latitude and height.

        Method to determine the pressure as a function of altitude and latitude,
        for calculating gaseous attenuation along an Earth-space path.
        This method is recommended when more reliable local data are not available.

        Parameters
        ----------
        lat : numpy.ndarray
            Latitude (degree)
        h : numpy.ndarray
            Height (km)
        season : str
            Season of the year (available values, 'summer', and 'winter').
            Default 'summer'


        Returns
        -------
        P: numpy.ndarray
            Pressure (hPa)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en
        """
        return self.instance.pressure(lat, h, season)

    def water_vapour_density(self, lat: float, h: np.ndarray, season: Optional[str] = 'summer') -> np.ndarray:
        """ Determine the water vapour density at a given latitude and height.

        Method to determine the water-vapour density as a
        function of altitude and latitude, for calculating gaseous attenuation
        along an Earth-space path. This method is recommended when more reliable
        local data are not available.


        Parameters
        ----------
        lat : numpy.ndarray
            Latitude (degree)
        h : numpy.ndarray
            Height (km)
        season : str
            Season of the year (available values, 'summer', and 'winter').
            Default 'summer'


        Returns
        -------
        rho: numpy.ndarray
            Water vapour density (g/m^3)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en
        """
        return self.instance.water_vapour_density(lat, h, season)

    def standard_temperature(self, h: np.ndarray, T_0: float) -> np.ndarray:
        """ Determine the standard temperature at a given height.

        Method to compute the temperature of an standard atmosphere at
        a given height. The reference standard atmosphere is based on the United
        States Standard Atmosphere, 1976, in which the atmosphere is divided into
        seven successive layers showing linear variation with temperature.


        Parameters
        ----------
        h : numpy.ndarray
            Height (km)
        T_0 : float
            Surface temperature (K)


        Returns
        -------
        T: numpy.ndarray
            Temperature (K)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en
        """
        return self.instance.standard_temperature(h, T_0)

    def standard_pressure(self, h: np.ndarray, T_0: float, P_0: float) -> np.ndarray:
        """ Determine the standard pressure at a given height.

        Method to compute the total atmopsheric pressure of an standard atmosphere
        at a given height.

        The reference standard atmosphere is based on the United States Standard
        Atmosphere, 1976, in which the atmosphere is divided into seven successive
        layers showing linear variation with temperature.


        Parameters
        ----------
        h : numpy.ndarray
            Height (km)
        T_0 : float
            Surface temperature (K)
        P_0 : float
            Surface pressure (hPa)


        Returns
        -------
        P: numpy.ndarray
            Pressure (hPa)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en
        """
        return self.instance.standard_pressure(h, T_0, P_0)

    def standard_water_vapour_density(self, h: np.ndarray, h_0: float, rho_0: float) -> np.ndarray:
        """ Determine the standard water vapour density at a given height.

        The reference standard atmosphere is based on the United States Standard
        Atmosphere, 1976, in which the atmosphere is divided into seven successive
        layers showing linear variation with temperature.


        Parameters
        ----------
        h : numpy.ndarray
            Height (km)
        h_0 : float
            Scale height (km)
        rho_0 : float
            Surface water vapour density (g/m^3)


        Returns
        -------
        rho: numpy.ndarray
            Water vapour density (g/m^3)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en
        """
        return self.instance.standard_water_vapour_density(h, h_0, rho_0)

    def standard_water_vapour_pressure(self, h: np.ndarray, h_0: Optional[float] = 2., rho_0: Optional[float] = 7.5) -> np.ndarray:
        """Determine the standard water vapour pressure at a given height.

        The reference standard atmosphere is based on the United States Standard
        Atmosphere, 1976, in which the atmosphere is divided into seven successive
        layers showing linear variation with temperature.


        Parameters
        ----------
        h : numpy.ndarray
            Height (km)
        h_0 : float
            Scale height (km)
            Default 2.
        rho_0 : float
            Surface water vapour density (g/m^3)
            Default 7.5


        Returns
        -------
        e: numpy.ndarray
            Water vapour pressure (hPa)


        References
        ----------
        [1] Reference Standard Atmospheres
        https://www.itu.int/rec/R-REC-P.835/en
        """

        return self.instance.standard_water_vapour_pressure(h, h_0, rho_0)

    def _get_season(self, lat: float, month: int) -> str:
        """_summary_

        Args:
            lat (float): _description_
            month (int): _description_

        Returns:
            str: _description_
        """
        if (lat > 20 and month in [6, 7, 8]) or (lat < -20 and month in [1, 2, 12]):
            season = 'summer'
        elif (lat > 20 and month in [1, 2, 12]) or (lat < -20 and month in [6, 7, 8]):
            season = 'winter'
        else:
            season = 'standard'

        return season

    def profile_extrapolation(self, lat: float, month: int, z: np.ndarray,
                              q: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Extrapolation of atmospheres to be used to determine 
        temperature, pressure and water-vapour pressure as a function 
        of altitude and latitude, for calculating gaseous attenuation when more reliable 
        local data are not available.

        Args:
            lat (float): Latitude (degree)
            month (int): Month of the year
            z (np.ndarray): Height (km)
            q (Tuple[np.ndarray, np.ndarray, np.ndarray]): Pressure (hPa), Temperature (K) and Relative humisity (frac) profiles

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Height, Pressure, Temperature, RH extrapolated profiles
        """
    
        if np.max(z) < 50:
            h_km = np.append(z, np.arange(max(z)+3.2, 50, 3.2))
            idx = np.where(self._height > np.max(h_km))
            self._height = np.append(h_km, self._height[idx])
        else:
            self._height = np.append(z, self._height[idx])

        season = self._get_season(lat, month)
        idx = np.where(self._height > np.max(z))

        p, t, rh = q

        if season != 'standard':
            tt = self.temperature(lat, self._height, season)
            pp = self.pressure(lat, self._height, season)
            wvd = self.water_vapour_density(lat, self._height, season)
        else:
            tt = self.instance.standard_temperature(self._height)
            pp = self.instance.standard_pressure(self._height)
            wvd = self.instance.standard_water_vapour_density(
                self._height)

        pres = np.append(p, pp[idx])
        temp = np.append(t, tt[idx])
        _rh = rho2rh(wvd, temp, pres)[0]
        rh = np.append(rh, _rh[idx])

        return self._height, pres, temp, rh


class _ITU835_6():
    """Class to model the ITU-R P.835_6 recommendation.

    The procedures to compute the reference standard atmosphere parameters
    pressented in these versions are identical to those included in version
    ITU_T P.835.

    Version:
       * P.835-6 (12/17) (Current version)
    """

    def __init__(self):
        self.__version__ = 6
        self.year = 2017
        self.month = 12
        self.link = 'https://www.itu.int/rec/R-REC-P.835-6-201712-I/en'
        self.h_km = np.hstack((np.linspace(0, 25, 26),  np.linspace(27.5, 50, 10),
                               np.linspace(55, 100, 10)))

    @staticmethod
    def standard_temperature(h, T_0=288.15):
        """

        """
        h_p = 6356.766 * h / (6356.766 + h)
        # Warnings because of sqrt are expected
        with np.errstate(invalid='ignore'):
            T = np.where(h_p <= 11, 288.15 - 6.5 * h_p,
                         np.where(np.logical_and(11 < h_p, h_p <= 20),
                                  216.65,
                                  np.where(np.logical_and(20 < h_p, h_p <= 32),
                                  216.65 + (h_p - 20),
                                  np.where(np.logical_and(32 < h_p, h_p <= 47),
                                           228.65 + 2.8 * (h_p - 32),
                                           np.where(np.logical_and(47 < h_p, h_p <= 51),
                                           270.65,
                                           np.where(np.logical_and(51 < h_p, h_p <= 71),
                                                    270.65 - 2.8 * (h_p - 51),
                                                    np.where(np.logical_and(71 < h_p, h_p <= 84.852),
                                                    214.65 - 2.0 * (h_p - 71),
                                                    np.where(np.logical_and(86 <= h, h <= 91),
                                                             186.8673,
                                                             np.where(np.logical_and(91 < h, h <= 100),
                                                             263.1905 - 76.3232 *
                                                             np.sqrt(
                                                                 (1 - ((h - 91)/19.9429)**2)),
                                                        195.08134)))))))))

        return T

    @staticmethod
    def standard_pressure(h, T_0=None, P_0=None):
        """

        """
        h_p = 6356.766 * h / (6356.766 + h)
        with np.errstate(invalid='ignore'):
            P = np.where(h_p <= 11,
                         1013.25 * (288.15 / (288.15 - 6.5 * h_p)
                                    )**(-34.1632 / 6.5),
                         np.where(np.logical_and(11 < h_p, h_p <= 20),
                                  226.3226 *
                                  np.exp(-34.1632 * (h_p - 11) / 216.65),
                                  np.where(np.logical_and(20 < h_p, h_p <= 32),
                                  54.74980 * (216.65 / (216.65 + (h_p - 20))
                                              ) ** 34.1632,
                                  np.where(np.logical_and(32 < h_p, h_p <= 47),
                                           8.680422 * (228.65 / (228.65 + 2.8 * (h_p - 32))) **
                                           (34.1632 / 2.8),
                                           np.where(np.logical_and(47 < h_p, h_p <= 51),
                                                    1.109106 *
                                                    np.exp(-34.1632 *
                                                           (h_p - 47) / 270.65),
                                                    np.where(np.logical_and(51 < h_p, h_p <= 71),
                                                             0.6694167 * (270.65 / (270.65 - 2.8 *
                                                                                    (h_p - 51)))**(-34.1632 / 2.8),
                                                             np.where(np.logical_and(71 < h_p, h_p <= 84.852),
                                                                      0.03956649 *
                                                                      (214.65 / (214.65 - 2.0 *
                                                                       (h_p - 71)))**(-34.1632 / 2.0),
                                                                      np.where(np.logical_and(86 <= h, h <= 100),
                                                                               np.exp(95.571899 - 4.011801 * h + 6.424731e-2 * h**2 -
                                                                                      4.789660e-4 * h**3 + 1.340543e-6 * h**4),
                                                                               1e-62)))))))).astype(float)

        return P

    @staticmethod
    def standard_water_vapour_density(h, h_0=2, rho_0=7.5):
        """

        """
        return rho_0 * np.exp(-h / float(h_0))

    def standard_water_vapour_pressure(self, h, h_0=2, rho_0=7.5):
        """

        """
        rho_h = self.standard_water_vapour_density(h, h_0, rho_0)
        T_h = self.standard_temperature(h)
        return rho_h * T_h / 216.7

    #  Low latitude standard atmosphere functions  (Section ITU-R P.835-6 2)  #
    @staticmethod
    def low_latitude_temperature(h):
        """Section 2 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h < 17)),
                        300.4222 - 6.3533 * h + 0.005886 * h**2,
                        np.where(np.logical_and((17 <= h), (h < 47)),
                        194 + (h - 17) * 2.533,
                        np.where(np.logical_and((47 <= h), (h < 52)), 270,
                                 np.where(np.logical_and((52 <= h), (h < 80)),
                                 270 - (h - 52) * 3.0714,
                            np.where(np.logical_and((80 <= h), (h <= 100)), 184, 184)))))

    def low_latitude_pressure(self, h):
        """Section 2 of Recommendation ITU-R P.835-6."""
        P10 = 284.8526   # Pressure at 10 km using the equation below
        P72 = 0.0313660  # Pressure at 72 km using the equation below
        return np.where(np.logical_and((0 <= h), (h <= 10)),
                        1012.0306 - 109.0338 * h + 3.6316 * h**2,
                        np.where(np.logical_and((10 < h), (h <= 72)),
                        P10 * np.exp(-0.147 * (h - 10)),
                        np.where(np.logical_and((72 < h), (h <= 100)),
                                 P72 * np.exp(-0.165 * (h - 72)), np.nan)))

    @staticmethod
    def low_latitude_water_vapour(h):
        """Section 3.1 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h <= 15)), 19.6542 *
                        np.exp(-0.2313 * h - 0.1122 * h**2 + 0.01351 * h**3 -
                               0.0005923 * h**4), 0)

    # Mid latitude standard atmosphere functions  (Section ITU-R P.835-6 3)
    @staticmethod
    def mid_latitude_temperature_summer(h):
        """Section 3.1 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h < 13)),
                        294.9838 - 5.2159 * h - 0.07109 * h**2,
                        np.where(np.logical_and((13 <= h), (h < 17)), 215.15,
                        np.where(np.logical_and((17 <= h), (h < 47)),
                                 215.15 * np.exp((h - 17) * 0.008128),
                                 np.where(np.logical_and((47 <= h), (h < 53)), 275,
                                 np.where(np.logical_and((53 <= h), (h < 80)),
                                          275 + 20 *
                                          (1 - np.exp((h - 53) * 0.06)),
                                          np.where(np.logical_and((80 <= h), (h <= 100)),
                                          175, np.nan))))))

    def mid_latitude_pressure_summer(self, h):
        """Section 3.1 of Recommendation ITU-R P.835-6."""
        P10 = 283.7096    # Pressure at 10 km using the equation below
        P72 = 0.03124022  # Pressure at 72 km using the equation below
        return np.where(
            np.logical_and((0 <= h), (h <= 10)),
            1012.8186 - 111.5569 * h + 3.8646 * h**2, np.where(
                np.logical_and((10 < h), (h <= 72)),
                P10 * np.exp(-0.147 * (h - 10)),
                np.where(
                    np.logical_and((72 < h), (h <= 100)),
                    P72 * np.exp(-0.165 * (h - 72)),
                    np.nan)))

    @staticmethod
    def mid_latitude_water_vapour_summer(h):
        """Section 3.1 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h <= 15)),
                        14.3542 * np.exp(-0.4174 * h - 0.02290 * h**2 +
                                         0.001007 * h**3), 0)

    @staticmethod
    def mid_latitude_temperature_winter(h):
        """Section 3.2 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h < 10)),
                        272.7241 - 3.6217 * h - 0.1759 * h**2,
                        np.where(np.logical_and((10 <= h), (h < 33)), 218,
                        np.where(np.logical_and((33 <= h), (h < 47)),
                                 218 + (h - 33) * 3.3571,
                                 np.where(np.logical_and((47 <= h), (h < 53)), 265,
                                 np.where(np.logical_and((53 <= h), (h < 80)),
                                          265 - (h - 53) * 2.0370,
                                          np.where(np.logical_and((80 <= h), (h <= 100)),
                                                   210, np.nan))))))

    def mid_latitude_pressure_winter(self, h):
        """Section 3.2 of Recommendation ITU-R P.835-6."""
        P10 = 258.9787    # Pressure at 10 km using the equation below
        P72 = 0.02851702  # Pressure at 72 km using the equation below
        return np.where(np.logical_and((0 <= h), (h <= 10)),
                        1018.8627 - 124.2954 * h + 4.8307 * h**2,
                        np.where(np.logical_and((10 < h), (h <= 72)),
                        P10 * np.exp(-0.147 * (h - 10)),
                        np.where(np.logical_and((72 < h), (h <= 100)),
                                 P72 * np.exp(-0.155 * (h - 72)), np.nan)))

    @staticmethod
    def mid_latitude_water_vapour_winter(h):
        """Section 3.2 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and(0 <= h, h <= 10),
                        3.4742 * np.exp(- 0.2697 * h - 0.03604 * h**2 +
                                        0.0004489 * h**3), 0)

    #  High latitude standard atmosphere functions  (Section ITU-R P.835-6 4)  #
    @staticmethod
    def high_latitude_temperature_summer(h):
        """Section 4.1 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h < 10)),
                        286.8374 - 4.7805 * h - 0.1402 * h**2,
                        np.where(np.logical_and((10 <= h), (h < 23)), 225,
                        np.where(np.logical_and((23 <= h), (h < 48)),
                                 225 * np.exp((h - 23) * 0.008317),
                                 np.where(np.logical_and((48 <= h), (h < 53)), 277,
                                 np.where(np.logical_and((53 <= h), (h < 79)),
                                          277 - (h - 53) * 4.0769,
                                          np.where(np.logical_and((79 <= h), (h <= 100)),
                                                   171, np.nan))))))

    def high_latitude_pressure_summer(self, h):
        """Section 4.1 of Recommendation ITU-R P.835-6."""
        P10 = 269.6138    # Pressure at 10 km using the equation below
        P72 = 0.04582115  # Pressure at 72 km using the equation below
        return np.where(np.logical_and((0 <= h), (h <= 10)),
                        1008.0278 - 113.2494 * h + 3.9408 * h**2,
                        np.where(np.logical_and((10 < h), (h <= 72)),
                        P10 * np.exp(-0.140 * (h - 10)),
                        np.where(np.logical_and((72 < h), (h <= 100)),
                                 P72 * np.exp(-0.165 * (h - 72)), np.nan)))

    @staticmethod
    def high_latitude_water_vapour_summer(h):
        """Section 4.1 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h <= 15)),
                        8.988 * np.exp(- 0.3614 * h - 0.005402 * h**2 -
                                       0.001955 * h**3), 0)

    @staticmethod
    def high_latitude_temperature_winter(h):
        """Section 4.2 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h < 8.5)),
                        257.4345 + 2.3474 * h - 1.5479 * h**2 +
                        0.08473 * h**3,
                        np.where(np.logical_and((8.5 <= h), (h < 30)), 217.5,
                        np.where(np.logical_and((30 <= h), (h < 50)),
                                 217.5 + (h - 30) * 2.125,
                                 np.where(np.logical_and((50 <= h), (h < 54)), 260,
                                 np.where(np.logical_and((54 <= h), (h <= 100)),
                                          260 - (h - 54) * 1.667, np.nan)))))

    def high_latitude_pressure_winter(self, h):
        """Section 4.2 of Recommendation ITU-R P.835-6."""
        P10 = 243.8718    # Pressure at 10 km using the equation below
        P72 = 0.02685355  # Pressure at 72 km using the equation below
        return np.where(np.logical_and((0 <= h), (h <= 10)),
                        1010.8828 - 122.2411 * h + 4.554 * h**2,
                        np.where(np.logical_and((10 < h), (h <= 72)),
                        P10 * np.exp(-0.147 * (h - 10)),
                        np.where(np.logical_and((72 < h), (h <= 100)),
                                 P72 * np.exp(-0.150 * (h - 72)), np.nan)))

    @staticmethod
    def high_latitude_water_vapour_winter(h):
        """Section 4.2 of Recommendation ITU-R P.835-6."""
        return np.where(np.logical_and((0 <= h), (h <= 10)),
                        1.2319 * np.exp(0.07481 * h - 0.0981 * h**2 +
                                        0.00281 * h**3), 0)

    def temperature(self, lat, h, season='summer'):
        """ Section 2 of Recommendation ITU-R P.835-6."""
        if season == 'summer':
            return np.where(
                np.abs(lat) < 22, self.low_latitude_temperature(h),
                np.where(
                    np.abs(lat) < 45, self.mid_latitude_temperature_summer(h),
                    self.high_latitude_temperature_summer(h)))
        elif season == 'winter':
            return np.where(
                np.abs(lat) < 22, self.low_latitude_temperature(h),
                np.where(
                    np.abs(lat) < 45, self.mid_latitude_temperature_winter(h),
                    self.high_latitude_temperature_winter(h)))
        else:
            raise ValueError("The value for argument 'season' is not correct."
                             "Valid values are 'summer' and 'winter'.")

    def pressure(self, lat, h, season='summer'):
        """ Section 2 of Recommendation ITU-R P.835-6."""
        if season == 'summer':
            return np.where(
                np.abs(lat) < 22, self.low_latitude_pressure(h),
                np.where(
                    np.abs(lat) < 45, self.mid_latitude_pressure_summer(h),
                    self.high_latitude_pressure_summer(h)))
        elif season == 'winter':
            return np.where(
                np.abs(lat) < 22, self.low_latitude_pressure(h),
                np.where(
                    np.abs(lat) < 45, self.mid_latitude_pressure_winter(h),
                    self.high_latitude_pressure_winter(h)))
        else:
            raise ValueError("The value for argument 'season' is not correct."
                             "Valid values are 'summer' and 'winter'")

    def water_vapour_density(self, lat, h, season='summer'):
        """ Section 2 of Recommendation ITU-R P.835-6."""
        if season == 'summer':
            return np.where(
                np.abs(lat) < 22, self.low_latitude_water_vapour(h),
                np.where(
                    np.abs(lat) < 45, self.mid_latitude_water_vapour_summer(h),
                    self.high_latitude_water_vapour_summer(h)))
        elif season == 'winter':
            return np.where(
                np.abs(lat) < 22, self.low_latitude_water_vapour(h),
                np.where(
                    np.abs(lat) < 45, self.mid_latitude_water_vapour_winter(h),
                    self.high_latitude_water_vapour_winter(h)))
        else:
            raise ValueError("The value for argument 'season' is not correct."
                             "Valid values are 'summer' and 'winter'")

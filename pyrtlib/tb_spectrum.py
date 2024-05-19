"""
Main class to compute Tb.
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

import warnings
from typing import Dict, Tuple, Union, Optional

import numpy as np
import pandas as pd

from .absorption_model import O2AbsModel, H2OAbsModel, N2AbsModel, LiqAbsModel
from .rt_equation import RTEquation


class TbCloudRTE(object):
    """
    Initialize TbCloudRTE
    """

    def __init__(self, z: np.ndarray, p: np.ndarray, t: np.ndarray, rh: np.ndarray, frq: np.ndarray,
                 angles: Optional[np.ndarray] = np.array([90.]),
                 o3n: Optional[np.ndarray] = None,
                 amu: Optional[Tuple] = None,
                 absmdl: Optional[str] = '',
                 ray_tracing: Optional[bool] = False,
                 from_sat: Optional[bool] = True,
                 cloudy: Optional[bool] = False):
        """Main class which computes brightness temperatures (Tb), mean
        radiating temperature (Tmr), and integrated absorption (Tau) for 
        clear or cloudy conditions. Also returns all integrated quantities
        that the original TBMODEL, Cyber Version, returned.  The input
        profiles are not modified within this subroutine.  It is assumed
        that the input profiles start at the antenna height (zX(1)).  The
        input profiles must reach 50.0 mb.  This subroutine uses the
        algorithms described in [Schroeder-Westwater-1991]_.

        Args:
            z (np.ndarray): Height profile (km).
            p (np.ndarray): Pressure profile (mb).
            t (np.ndarray): Temperature profile (K).
            rh (np.ndarray): Relative humidity profile (fraction).
            frq (np.ndarray): Channel frequencies (GHz).
            angles (Optional[np.ndarray], optional): Elevation anglesX (deg).. Defaults to 90.
            o3n (Optional[np.ndarray], optional): Ozone profile. Defaults to None.
            amu (Optional[Tuple], optional): Absorption model uncertainties. Defined by :py:func:`~pyrtlib.uncertainty.SpectroscopicParameter`. Defaults to None.
            absmdl (Optional[str], optional): Absorption model. Defaults to ''.
            ray_tracing (Optional[bool], optional): Wether True it computes ray tracing for
                                            distance between layers, otherwise use simple plane
                                            parallel assumption. Defaults to False.
            from_sat (Optional[bool], optional): Wether True (default) compute upwelling Tb, 
                                            otherwise downwelling Tb are computed. Defaults to True.
            cloudy (Optional[bool], optional): Wether True CLW must be passed. Defaults to False.
        """

        self.z = z
        self.p = p
        self.tk = t
        self.rh = rh
        self.frq = frq
        self.angles = angles
        self.o3n = o3n
        self.amu = amu

        dz = np.ediff1d(self.z)
        if all(dz) > 0:
            pass
        elif all(dz) < 0:
            self.z = z[::-1]
            self.p = p[::-1]
            self.tk = t[::-1]
            self.rh = rh[::-1]
            print("Profile flipped up/douwn")
        else:
            raise SystemExit("ERROR: input profile seems incorrect. "
                             "It must be monotonically increasing or decreasing")

        if len(self.p) < 25 or min(self.p) >= 10:
            warnings.warn(f"Number of levels too low ({len(self.p)}) or "
                          f"minimum pressure value lower than 10 hPa ({min(self.p)}). "
                          "Please considering profile extrapolation. Levels number must be higher than 25 "
                          "and pressure value lower than 10 hPa")

        self.ray_tracing = ray_tracing
        self._satellite = from_sat
        self.cloudy = cloudy
        self._uncertainty = True if self.amu else False

        self.nl = len(z)
        self.nf = len(frq)
        self.nang = len(angles)

        self.ice = False
        # set emissivity
        self._es = np.repeat(1.0, self.nf)

        # convert height profile to (km above antenna height)
        self.z0 = self.z[0]
        self.z -= self.z0

        # Allocation
        self.sptaudry = np.zeros((self.nf, self.nang))
        self.sptauwet = np.zeros((self.nf, self.nang))
        self.sptauliq = np.zeros((self.nf, self.nang))
        self.sptauice = np.zeros((self.nf, self.nang))
        self.awet = np.zeros((self.nf, self.nang, self.nl))
        self.adry = np.zeros((self.nf, self.nang, self.nl))
        self.ptaudry = np.zeros((self.nf, self.nang, self.nl))
        self.ptaulay = np.zeros((self.nf, self.nang, self.nl))
        self.ptauwet = np.zeros((self.nf, self.nang, self.nl))
        self.ptauliq = np.zeros((self.nf, self.nang, self.nl))
        self.ptauice = np.zeros((self.nf, self.nang, self.nl))
        self.tbtotal = np.zeros((self.nf, self.nang))
        self.tbatm = np.zeros((self.nf, self.nang))
        self.tmr = np.zeros((self.nf, self.nang))
        self.tmrcld = np.zeros((self.nf, self.nang))
        self.srho = np.zeros((1, self.nang))
        self.swet = np.zeros((1, self.nang))
        self.sdry = np.zeros((1, self.nang))
        self.sliq = np.zeros((1, self.nang))
        self.sice = np.zeros((1, self.nang))

        if absmdl:
            self.init_absmdll(absmdl)

    @property
    def satellite(self) -> bool:
        """If :code:`True` computes an upward-propagating brightness-temperature spectrum
        otherwise a downward-propagating brightness-temperature 
        spectrum at the bottom of the atmosphere will be performed.
        """
        return self._satellite

    @satellite.setter
    def satellite(self, sat: bool) -> None:
        if isinstance(sat, bool):
            self._satellite = sat
        else:
            raise ValueError("Please enter True or False")

    @property
    def emissivity(self) -> Union[float, np.ndarray]:
        """Surface emissivity. Default to 1.0
        """
        return self._es

    @emissivity.setter
    def emissivity(self, emissivity: Union[float, np.ndarray]) -> None:
        """Setter for surface emissivty

        Args:
            emissivity (np.float): [description]

        Raises:
            ValueError: [description]
        """
        if isinstance(emissivity, float):
            self._es = np.repeat(emissivity, self.nf)
        elif isinstance(emissivity, np.ndarray):
            self._es = emissivity
        else:
            raise ValueError(
                "Please enter a valid value or array for emissivity")

    def set_amu(self, amu: Tuple) -> None:
        """Set absorption model uncertainties

        Args:
            amu (Dict): The spectroscopic parameters dictionary

        See also:
            :py:func:`~pyrtlib.uncertainty.SpectroscopicParameter`
        """
        self.amu = amu
        self._uncertainty = True

    def init_absmdl(self, absmdl: str):
        """Initialize absorption models.

        Args:
            absmdl (str): Absorption model.
        """
        # Defines models
        H2OAbsModel.model = absmdl
        O2AbsModel.model = absmdl
        N2AbsModel.model = absmdl
        LiqAbsModel.model = absmdl

    def _init_linelist(self):
        H2OAbsModel.set_ll()
        O2AbsModel.set_ll()

    def init_cloudy(self, cldh: np.ndarray, denice: np.ndarray, denliq: np.ndarray) -> None:
        """Initialize cloudy conditions parameters.

        Args:
            cldh (numpy.ndarray): Cloud base/top heights (km MSL)
            denice (numpy.ndarray): Ice density profile (:math:`g/m^3`).
            denliq (numpy.ndarray): Liquid density profile (:math:`g/m^3`).
        """
        # convert cloud base and cloud top to (km above antenna height)
        # compute (beglev) and (endlev)
        ncld = cldh.shape[1]
        self.cldh = cldh - self.z0
        self.beglev = np.zeros((ncld,))
        self.endlev = np.zeros((ncld,))
        if getattr(self, 'cloudy'):
            for j in range(0, ncld):
                for i in range(0, self.nl):
                    if self.z[i] == self.cldh[0, j]:
                        self.beglev[j] = i
                    if self.z[i] == self.cldh[1, j]:
                        self.endlev[j] = i

            self.denice = denice
            self.denliq = denliq
        else:
            warnings.warn("It seems that TbCloudRTE.cloudy attribute is not set to True. "
                          "Sets it to True for running model in cloudy condition.")
            # raise AttributeError("Set cloudy to True before running init_cloudy()")

    def execute(self, only_bt: bool = True) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, np.ndarray]]]:
        """This function computes Brightness Temperature and other radiometric parameters.

        Args:
            only_bt (bool): If True returns only brightness temperature. Default to True.

        Returns a pandas dataframe containing:

            tbtotal:
                brightness temperature (K) includes cosmic background; 
                indexed by frequency and elevation angle

            tbatm:
                atmospheric brightness temperature (K), no cosmic; 
                background;indexed by frequency and elevation angle

            tmr:
                mean radiating temperature of the atmosphere (K);
                indexed by frequency and elevation angle

            tmrcld:
                mean radiating temperature (K) of the lowest cloud layer;
                indexed by frequency and elevation angle

            taudry:
                dry air absorption integrated over each ray path (Np);
                indexed by frequency and elevation angle

            tauwet:
                water vapor absorption integrated over each ray path (Np);
                indexed by frequency and elevation angle

            tauliq:
                cloud liquid absorption integrated over each ray path (Np);
                indexed by frequency and elevation angle

            tauice:
                cloud ice absorption integrated over each ray path (Np);
                indexed by frequency and elevation angle

        and with `only_bt=False` also a dictionary containing:

            taulaydry:
                dry air absorption integrated over each ray path (Np);
                indexed by frequency, elevation angle and height profile

            taulaywet:
                water vapor absorption integrated over each ray path (Np);
                indexed by frequency, elevation angle and height profile

            taulayliq:
                cloud liquid absorption integrated over each ray path (Np);
                indexed by frequency, elevation angle and height profile

            taulayice:
                cloud ice absorption integrated over each ray path (Np);
                indexed by frequency, elevation angle and height profile

            srho:
                water vapor density integrated along each ray path (cm);
                indexed by elevation angle

            swet:
                wet refractivity integrated along each ray path (cm);
                indexed by elevation angle

            sdry:
                dry refractivity integrated along each ray path (cm);
                indexed by elevation angle

            sliq:
                cloud ice density integrated along each ray path (cm);
                indexed by elevation angle

            sice:
                cloud liquid density integrated along each ray path (cm);
                indexed by elevation angle

        Returns:
            Union[pandas.DataFrame, Tuple[pandas.DataFrame, Dict[str, numpy.ndarray]]]: A pandas Dataframe with the brigthness temperature simulated.
            If only_bt = False it also returns all intermediate RT variables.
        """

        self._init_linelist()

        if self.cloudy and LiqAbsModel.model in ['R98', 'R03']:
            warnings.warn(
                "Model {} for liquid cloud absorption is outdated. "
                "Use the more recent model by Rosenkranz, 2015, from R16 model onwards".format(LiqAbsModel.model), stacklevel=2)

        # Set RTE
        RTEquation._from_sat = self._satellite

        # compute vapor pressure and vapor density
        e, rho = RTEquation.vapor(self.tk, self.rh, self.ice)
        # compute refractivity
        dryn, wetn, refindx = RTEquation.refractivity(self.p, self.tk, e)
        for k in range(0, self.nang):
            # Compute distance between each level (ds)
            if self.ray_tracing:
                ds = RTEquation.ray_tracing(
                    self.z, refindx, self.angles[k], self.z0)
            else:
                amass = 1 / np.sin(self.angles[k] * np.pi / 180)
                ds = np.append([0], np.diff(self.z) * amass)
            # ds = [0; diff(z)]; # in alternative simple diff of z

            # Integrate over path (ds)
            self.srho[:, k], _ = RTEquation.exponential_integration(
                True, rho, ds, 0, self.nl, 0.1)
            self.swet[:, k], _ = RTEquation.exponential_integration(
                True, wetn, ds, 0, self.nl, 0.1)
            self.sdry[:, k], _ = RTEquation.exponential_integration(
                True, dryn, ds, 0, self.nl, 0.1)
            if self.cloudy:
                self.sliq[:, k] = RTEquation.cloud_integrated_density(
                    self.denliq, ds, self.beglev, self.endlev)
                self.sice[:, k] = RTEquation.cloud_integrated_density(
                    self.denice, ds, self.beglev, self.endlev)

            # handle each frequency
            # this are based on NOAA RTE fortran routines
            for j in range(0, self.nf):
                RTEquation._emissivity = self._es[j]
                # Rosenkranz, personal communication, 2019/02/12 (email)
                self.awet[j, k, :], self.adry[j, k, :] = RTEquation.clearsky_absorption(self.p, self.tk, e, self.frq[j],
                                                            self.o3n, self.amu if self._uncertainty else None)
                self.sptauwet[j, k], \
                    self.ptauwet[j, k, :] = RTEquation.exponential_integration(
                        True, self.awet[j, k, :], ds, 1, self.nl, 1)
                self.sptaudry[j, k], \
                    self.ptaudry[j, k, :] = RTEquation.exponential_integration(
                        True, self.adry[j, k, :], ds, 1, self.nl, 1)
                if self.cloudy:
                    aliq, aice = RTEquation.cloudy_absorption(
                        self.tk, self.denliq, self.denice, self.frq[j])
                    self.sptauliq[j, k], \
                        self.ptauliq[j, k, :] = RTEquation.exponential_integration(
                            False, aliq, ds, 1, self.nl, 1)
                    self.sptauice[j, k], \
                        self.ptauice[j, k, :] = RTEquation.exponential_integration(
                            False, aice, ds, 1, self.nl, 1)

                self.ptaulay[j, k, :] = self.ptauwet[j, k, :] + \
                    self.ptaudry[j, k, :] + \
                    self.ptauice[j, k, :] + \
                    self.ptauliq[j, k, :]
                # [boftotl,boftatm,boftmr,PSPtauprof,hvk] = Planck_xxx(frq(j),tk,Ptaulay(j,k,:));
                boftotl, boftatm, boftmr, psp_tauprof, hvk, _, _ = RTEquation.planck(self.frq[j], self.tk,
                                                                                     self.ptaulay[j, k, :])
                if self.cloudy:
                    self.tmrcld[j, k] = RTEquation.cloud_radiating_temperature(self.beglev[0], self.endlev[0], hvk,
                                                                               psp_tauprof,
                                                                               boftatm)
                # assign output values
                self.tbtotal[j, k] = RTEquation.bright(hvk, boftotl)
                self.tbatm[j, k] = RTEquation.bright(hvk, boftatm[self.nl - 1])
                self.tmr[j, k] = RTEquation.bright(hvk, boftmr)

        df = pd.DataFrame()
        for i, a in enumerate(self.angles):
            df_i = pd.DataFrame({'tbtotal': self.tbtotal.T[i],
                                 'tbatm': self.tbatm.T[i],
                                 'tmr': self.tmr.T[i],
                                 'tmrcld': self.tmrcld.T[i],
                                 'tauwet': self.sptauwet.T[i],
                                 'taudry': self.sptaudry.T[i],
                                 'tauliq': self.sptauliq.T[i],
                                 'tauice': self.sptauice.T[i],
                                 'angle': np.full((len(self.frq),), a)
                                 })
            df = pd.concat([df, df_i])

        if only_bt:
            return df
        else:
            return df, {'taulaywet': self.ptauwet, 'taulaydry': self.ptaudry,
                        'taulayliq': self.ptauliq, 'taulayice': self.ptauice,
                        'awet': self.awet, 'adry': self.adry,
                        'srho': self.srho, 'swet': self.swet, 'sdry': self.sdry,
                        'sliq': self.sliq, 'sice': self.sice}

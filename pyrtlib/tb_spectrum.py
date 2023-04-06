"""
Main script !!!(right now only for testing)!!!.
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
from .utils import import_lineshape


class TbCloudRTE(object):
    """
    Initialize TbCloudRTE
    """

    def __init__(self, z: np.ndarray, p: np.ndarray, tk: np.ndarray, rh: np.ndarray, frq: np.ndarray,
                 angles: Optional[np.ndarray] = np.array([90.]),
                 o3n: Optional[np.ndarray] = None,
                 amu: Optional[Tuple] = None,
                 absmdl: Optional[str] = '',
                 ray_tracing: Optional[bool] = True,
                 from_sat: Optional[bool] = True,
                 cloudy: Optional[bool] = False):
        """User interface which computes brightness temperatures (Tb), mean
        radiating temperature (Tmr), and integrated absorption (Tau) for 
        clear or cloudy conditions.  Also returns all integrated quantities
        that the original TBMODEL, Cyber Version, returned.  The input
        profiles are not modified within this subroutine.  It is assumed
        that the input profiles start at the antenna height (zX(1)).  The
        input profiles must reach 50.0 mb.  This subroutine uses the
        algorithms described in Schroeder and Westwater (1991).

        Args:
            z (numpy.ndarray): height profile (km MSL).
            p (numpy.ndarray): pressure profile (mb).
            tk (numpy.ndarray): temperature profile (K).
            rh (numpy.ndarray): relative humidity profile (fraction).
            frq (numpy.ndarray): channel frequencies (GHz).
            angles (numpy.ndarray): elevation anglesX (deg).
            absmdl (str, optional): absorption model for WV. Defaults to ''.
            ray_tracing (bool, optional):   if True (default) it computes ray tracing (RayTrac_xxx) for
                                            distance between layers; otherwise use simple plane
                                            parallel assumption (i.e. ds = diff(z)*airmass;. Defaults to True.
            from_sat (bool, optional): [description]. Defaults to True.
            cloudy (bool, optional): [description]. Defaults to True.
        """

        self.z = z
        self.p = p
        self.tk = tk
        self.rh = rh
        self.frq = frq
        self.angles = angles
        self.o3n = o3n
        self.amu = amu

        self.ray_tracing = ray_tracing
        self._satellite = from_sat
        self.cloudy = cloudy
        self._uncertainty = False

        self.nl = len(z)
        self.nf = len(frq)
        self.nang = len(angles)

        self.ice = False
        # set emissivity
        self._es = np.repeat(1.0, self.nf)

        # ... convert height profile to (km above antenna height) ...
        self.z0 = self.z[0]
        self.z -= self.z0

        # Allocation
        self.sptaudry = np.zeros((self.nf, self.nang))
        self.sptauwet = np.zeros((self.nf, self.nang))
        self.sptauliq = np.zeros((self.nf, self.nang))
        self.sptauice = np.zeros((self.nf, self.nang))
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
        """If :code:`True` performing model calculation from satellite otherwise from ground

        Returns:
            bool: If True performing calculation from satellite 
                    otherwise from ground. Default to True.
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

        Returns:
            np.float: The surface emissivity.
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
            raise ValueError("Please enter a valid value or array for emissivity")

    def set_amu(self, amu: Tuple) -> Tuple:
        self.amu = amu

    def init_absmdl(self, absmdl: str):
        """Initialize absorption models.

        Args:
            absmdl (str): Absorption model for WV
        """
        if absmdl == 'uncertainty':
            O2AbsModel.model = 'rose18'
            O2AbsModel.o2ll = import_lineshape('o2ll_{}'.format('rose18'))
            H2OAbsModel.model = 'rose17'
            H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format('rose17'))
            N2AbsModel.model = 'rose03'
            LiqAbsModel.model = 'rose16'
            self._uncertainty = True
        else:
            # Defines models
            try:
                O2AbsModel.model = absmdl
                O2AbsModel.o2ll = import_lineshape('o2ll_{}'.format(absmdl))
            except KeyError as e:
                warnings.warn("The lines list {} was not found. You have to define absorption model manually".format(e))
            try:
                H2OAbsModel.model = absmdl
                H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(absmdl))
            except KeyError as e:
                warnings.warn("The lines list {} was not found".format(e))
            
            N2AbsModel.model = absmdl
            LiqAbsModel.model = absmdl

    def init_cloudy(self, cldh: np.ndarray, denice: np.ndarray, denliq: np.ndarray) -> None:
        """Initialize cloudy conditions parameters.

        Args:
            cldh (numpy.ndarray): cloud base/top heights (km MSL
            denice (numpy.ndarray): ice density profile (g/m**3); density fraction = 1.0
            denliq (numpy.ndarray): liquid density profile (g/m**3); density fraction = 1.0
        """
        # ... convert cloud base and cloud top to (km above antenna height) ...
        # ... compute (beglev) and (endlev) ...
        ncld = cldh.shape[1]
        self.cldh = cldh - self.z0
        self.beglev = np.zeros((ncld,))
        self.endlev = np.zeros((ncld,))
        if getattr(self, 'cloudy'):
            for j in range(0, ncld):
                for i in range(0, self.nl):
                    if self.z[i] == self.cldh[0, j]: self.beglev[j] = i
                    if self.z[i] == self.cldh[1, j]: self.endlev[j] = i

            self.denice = denice
            self.denliq = denliq
        else:
            warnings.warn("It seems that TbCloudRTE.cloudy attribute is not set to True. "
                          "Sets it to True for running model in cloudy condition.")
            # raise AttributeError("Set cloudy to True before running init_cloudy()")

    def execute(self, only_bt: bool = True) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, np.ndarray]]]:
        """Execution of main script.

        Args:
            only_bt (bool): If True (default) returns only brightness temperature. Default to True.
        
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

            and a dictionary containing:

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
            Union[pandas.DataFrame, Tuple[pandas.DataFrame, Dict[str, numpy.ndarray]]]: [description]
        """

        # Set RTE
        RTEquation.from_sat = self._satellite

        # ... compute vapor pressure and vapor density ...
        e, rho = RTEquation.vapor(self.tk, self.rh, self.ice)
        # ... compute refractivity ...
        dryn, wetn, refindx = RTEquation.refractivity(self.p, self.tk, e)
        for k in range(0, self.nang):
            # ... Compute distance between each level (ds) ...
            if self.ray_tracing:
                ds = RTEquation.ray_tracing(self.z, refindx, self.angles[k], self.z0)
            else:
                amass = 1 / np.sin(self.angles[k] * np.pi / 180)
                ds = np.append([0], np.diff(self.z) * amass)
            # ds = [0; diff(z)]; # in alternative simple diff of z

            # ... Integrate over path (ds) ...
            self.srho[k], _ = RTEquation.exponential_integration(True, rho, ds, 0, self.nl, 0.1)
            self.swet[k], _ = RTEquation.exponential_integration(True, wetn, ds, 0, self.nl, 0.1)
            self.sdry[k], _ = RTEquation.exponential_integration(True, dryn, ds, 0, self.nl, 0.1)
            if self.cloudy:
                self.sliq[k] = RTEquation.cloud_integrated_density(self.denliq, ds, self.beglev, self.endlev)
                self.sice[k] = RTEquation.cloud_integrated_density(self.denice, ds, self.beglev, self.endlev)

            # ... handle each frequency ...
            # this are based on NOAA RTE fortran routines
            for j in range(0, self.nf):
                RTEquation.emissivity = self._es[j]
                if self._uncertainty:
                    awet, adry = RTEquation.clearsky_absorption_uncertainty(self.p, self.tk, e, self.frq[j], self.amu)
                else:
                    # Rosenkranz, personal communication, 2019/02/12 (email)
                    awet, adry = RTEquation.clearsky_absorption(self.p, self.tk, e, self.frq[j], self.o3n)
                self.sptauwet[j, k], \
                self.ptauwet[j, k, :] = RTEquation.exponential_integration(1, awet, ds, 0, self.nl, 1)
                self.sptaudry[j, k], \
                self.ptaudry[j, k, :] = RTEquation.exponential_integration(1, adry, ds, 0, self.nl, 1)
                if self.cloudy:
                    aliq, aice = RTEquation.cloudy_absorption(self.tk, self.denliq, self.denice, self.frq[j])
                    self.sptauliq[j, k], \
                    self.ptauliq[j, k, :] = RTEquation.exponential_integration(0, aliq, ds, 0, self.nl, 1)
                    self.sptauice[j, k], \
                    self.ptauice[j, k, :] = RTEquation.exponential_integration(0, aice, ds, 0, self.nl, 1)

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
                # ... assign output values ...
                self.tbtotal[j, k] = RTEquation.bright(hvk, boftotl)
                self.tbatm[j, k] = RTEquation.bright(hvk, boftatm[self.nl - 1])
                self.tmr[j, k] = RTEquation.bright(hvk, boftmr)

        df = pd.DataFrame({'tbtotal': self.tbtotal.T[0],
                           'tbatm': self.tbatm.T[0],
                           'tmr': self.tmr.T[0],
                           'tmrcld': self.tmrcld.T[0],
                           'tauwet': self.sptauwet.T[0],
                           'taudry': self.sptaudry.T[0],
                           'tauliq': self.sptauliq.T[0],
                           'tauice': self.sptauice.T[0]
                           })

        if only_bt:
            return df
        else:
            return df, {'taulaywet': self.ptauwet, 'taulaydry': self.ptaudry,
                        'taulayliq': self.ptauliq, 'taulayice': self.ptauice,
                        'srho': self.srho, 'swet': self.swet, 'sdry': self.sdry,
                        'sliq': self.sliq, 'sice': self.sice}

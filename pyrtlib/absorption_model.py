# -*- coding: utf-8 -*-
"""
This class contains the absorption model used in pyrtlib.
"""

import types
import os
from typing import Tuple, Union, List, Optional, Dict
from netCDF4 import Dataset

import numpy as np

from pyrtlib.climatology import AtmosphericProfiles as atmp
from .utils import dilec12, _dcerror, constants, gas_mass, import_lineshape

PATH = os.path.dirname(os.path.abspath(__file__))


class AbsModelError(Exception):
    """Exception raised for errors in the input model.

    Attributes:
        model -- input model which caused the error
        message -- explanation of the error
    """

    def __init__(self, model, message):
        self.model = model
        self.message = message

        super().__init__(self.message)


class AbsModelError(Exception):
    """Exception raised for errors in the input model.

    Attributes:
        model -- input model which caused the error
        message -- explanation of the error
    """

    def __init__(self, model, message):
        self.model = model
        self.message = message

        super().__init__(self.message)


class AbsModel:
    """This is an abstraction class to define the absorption model.
    """

    def __init__(self) -> None:
        super(AbsModel, self).__init__()
        self._model = ''

    @property
    def model(self) -> str:
        """Getter/Setter for absorption model"""
        return self._model

    @model.setter
    def model(self, model: str) -> None:
        """Setter for absorption model"""
        if model and isinstance(model, str):
            self._model = model
        else:
            raise ValueError("Please enter a valid absorption model")

    @staticmethod
    def set_ll() -> None:
        """Set the linelist to the absorption model.

        See also:
            :py:func:`~pyrtlib.utils.import_lineshape`

        Example:
            .. code-block:: python

                from pyrtlib.absorption_model import H2OAbsModel
                H2OAbsModel.model = 'R21SD'
                H2OAbsModel.set_ll()

        .. note::
            Model must be set with `model` property before calling this method (see Example).
        """

    @staticmethod
    def implemented_models() -> Dict[str, List[str]]:
        """Return all the implemented absorption models.

        Returns:
            List[str]: The list the implemented absorption models

        Example:
            .. code-block:: python

                >>> from pyrtlib.absorption_model import AbsModel
                >>> AbsModel.implemented_models()
                {'Oxygen': ['R98',
                    'R03',
                    'R16',
                    'R17',
                    'R18',
                    'R19',
                    'R19SD',
                    'R20',
                    'R20SD',
                    'R22'],
                    'WaterVapour': ['R98',
                    'R03',
                    'R16',
                    'R17',
                    'R18',
                    'R19',
                    'R19SD',
                    'R20',
                    'R20SD',
                    'R21SD',
                    'R22SD'],
                    'Ozone': ['R18', 'R22']}
        """
        oxygen_netcdf = Dataset(os.path.join(
            PATH, '_lineshape', 'o2_lineshape.nc'), mode='r')
        wv_netcdf = Dataset(os.path.join(
            PATH, '_lineshape', 'h2o_lineshape.nc'), mode='r')
        ozone_netcdf = Dataset(os.path.join(
            PATH, '_lineshape', 'o3_lineshape.nc'), mode='r')

        model = {"Oxygen": list(oxygen_netcdf.groups.keys()),
                 "WaterVapour": list(wv_netcdf.groups.keys()),
                 "Ozone": list(ozone_netcdf.groups.keys())}

        return model


class LiqAbsModel(AbsModel):
    """This class contains the absorption model used in pyrtlib.
    """

    @staticmethod
    def liquid_water_absorption(water: np.ndarray, freq: np.ndarray, temp: np.ndarray) -> np.ndarray:
        """Computes absorption in Nepers/km by suspended liquid water droplets.

        Args:
            water (numpy.ndarray): Liquid water content (:math:`g/m^3`) - (mass of liquid water per volume of dry air).
            freq (numpy.ndarray): Frequency (GHz) - (valid from 0 to 1000 GHz).
            temp (numpy.ndarray): Temperature (K).

        Returns:
            np.ndarray: Liquid water particels absorption (Np/km)

        References
        ----------
        .. [1] [Liebe-Hufford-Manabe]_.
        .. [2] [Liebe-Hufford-Cotton]_.
        .. [3] [Rosenkranz-1988]_.
        .. [4] [Rosenkranz-2015]_.

        .. note::
            Revision history:

            * PWR 08/03/92 Original Version
            * PWR 12/14/98 Temp dependence of eps2 eliminated to agree with MPM93 
            * PWR 06/05/15 Using dilec12 for complex dielectric constant
        """
        if water <= 0:
            abliq = 0
            return abliq

        if LiqAbsModel.model in ['R03', 'R98']:
            theta1 = 1.0 - 300.0 / temp
            eps0 = 77.66 - np.dot(103.3, theta1)
            eps1 = np.dot(0.0671, eps0)
            eps2 = 3.52
            fp = np.dot((np.dot(316.0, theta1) + 146.4), theta1) + 20.2
            if LiqAbsModel.model == 'R03':
                fp = 20.1 * np.exp(7.88 * theta1)
            fs = np.dot(39.8, fp)
            eps = (eps0 - eps1) / complex(1.0, freq / fp) + \
                (eps1 - eps2) / complex(1.0, freq / fs) + eps2
        elif LiqAbsModel.model in ['R17', 'R16', 'R19', 'R20', 'R19SD', 'R22SD', 'R23', 'R23SD', 'R24']:
            eps = dilec12(freq, temp)
        else:
            raise ValueError(
                '[AbsLiq] No model available with this name: {} . Sorry...'.format(LiqAbsModel.model))

        re = (eps - 1.0) / (eps + 2.0)

        abliq = -0.06286 * np.imag(re) * freq * water

        return abliq


class N2AbsModel(AbsModel):
    """This class contains the absorption model used in pyrtlib.
    """

    @staticmethod
    def n2_absorption(t: np.ndarray, p: np.ndarray, f: np.ndarray) -> np.ndarray:
        """Collision-Induced Power Absorption Coefficient (Neper/km) in air
        with modification of 1.34 to account for O2-O2 and O2-N2 collisions, as calculated by [Boissoles-2003]_.

        Args:
            t (numpy.ndarray): Temperature (K).
            p (numpy.ndarray): Dry air pressure (mb).
            f (numpy.ndarray): Frequency (GHz) - (valid 0-2000 GHz).

        Raises:
            ValueError: Raises to error whether inputn model is unorrect or not available

        Returns:
            np.ndarray: Nitrogen continum absorption terms in Np/km

        References
        ----------
        .. [1] See eq. 2.6 in [Mätzler-Rosenkranz-2006]_.
        .. [2] [Borysow-Frommhold-1986]_.
        .. [3] [Boissoles-2003]_.
        """

        th = 300.0 / t
        fdepen = 0.5 + 0.5 / (1.0 + (f / 450.0) ** 2)
        if N2AbsModel.model in ['R16', 'R17', 'R18', 'R19', 'R19SD']:
            l, m, n = 6.5e-14, 3.6, 1.34
        elif N2AbsModel.model in ['R20', 'R20SD', 'R21SD', 'R22', 'R22SD', 'R23', 'R23SD', 'R24']:
            l, m, n = 9.95e-14, 3.22, 1
        elif N2AbsModel.model == 'R03':
            l, m, n = 6.5e-14, 3.6, 1.29
        elif N2AbsModel.model in ['R98']:
            l, m, n = 6.4e-14, 3.55, 1
            fdepen = 1
        else:
            raise ValueError(
                '[AbsN2] No model available with this name: {} . Sorry...'.format(N2AbsModel.model))

        bf = l * fdepen * p * p * f * f * th ** m

        abs_n2 = n * bf

        return abs_n2


class H2OAbsModel(AbsModel):
    """This class contains the :math:`H_2O` absorption model used in pyrtlib.
    """

    def __init__(self) -> None:
        super(H2OAbsModel, self).__init__()
        self._h2oll = None

    @property
    def h2oll(self) -> types.MethodType:
        """Getter for line list coefficient"""
        return self._h2oll

    @h2oll.setter
    def h2oll(self, h2oll) -> None:
        """Setter for line list coefficient"""
        if h2oll and isinstance(h2oll, types.MethodType):
            self._h2oll = h2oll
        else:
            raise ValueError("Please enter a valid absorption model")

    @staticmethod
    def set_ll() -> None:
        if H2OAbsModel.model not in H2OAbsModel.implemented_models()['WaterVapour']:
            raise AbsModelError(
                H2OAbsModel.model, f"Model {H2OAbsModel.model} is not available. It is necessary to define water vapour absorption model manually")
        H2OAbsModel.h2oll = import_lineshape("h2oll")

    def h2o_continuum(self, frq: np.ndarray, vx: np.ndarray, nfreq: int):
        nf = 6
        deltaf = 299.792458
        selfcon = np.array([2.877e-21, 2.855e-21, 2.731e-21,
                           2.49e-21, 2.178e-21, 1.863e-21])
        selftexp = np.array([6.413, 6.414, 6.275, 6.049, 5.789, 5.557])
        t = 300.0 / vx
        theta = 296./t

        a = np.zeros((nf))
        cs = np.zeros((nfreq))
        for j in range(nf):
            a[j] = 6.532e+12*selfcon[j]*theta**(selftexp[j]+3.)
        a = np.insert(a, 0, a[1], axis=0)

        for i in range(nfreq):
            fj = frq/deltaf
            j = int(fj)
            j = np.minimum(j, nf-2)
            p = fj - j
            c = (3.-2.*p)*p*p
            b = 0.5*p*(1.-p)
            b1 = b*(1.-p)
            b2 = b*p
            cs[i] = -a[j]*b1+a[j+1]*(1.-c+b2)+a[j+2]*(c+b1)-a[j+3]*b2

        return cs

    def h2o_absorption(self, pdrykpa: np.ndarray, vx: np.ndarray, ekpa: np.ndarray, frq: np.ndarray, amu: Optional[dict] = None) -> Union[
            Tuple[np.ndarray, np.ndarray], None]:
        """Compute absorption coefficients in atmosphere due to water vapor for all models.

        Args:
            pdrykpa (numpy.ndarray): Dry air pressure (kPa).
            vx (numpy.ndarray): Theta (adim) - (normalised temperature 300/t(K)).
            ekpa (numpy.ndarray): Water vapor partial pressure (kPa).
            frq (numpy.ndarray): Frequency (GHz) - (valid 0-1000 GHz).

        Returns:
            Union[ Tuple[numpy.ndarray, numpy.ndarray], None]: WV line and continuum absorption terms (ppm)

        References
        ----------
        .. [1] [Rosenkranz-2017]_.

        Example:
            .. code-block:: python

                import numpy as np
                from pyrtlib.rt_equation import RTEquation
                from pyrtlib.absorption_model import H2OAbsModel, AbsModel
                from pyrtlib.climatology import AtmosphericProfiles as atmp
                from pyrtlib.utils import ppmv2gkg, mr2rh, import_lineshape

                z, p, d, tk, md = atmp.gl_atm(atmp.TROPICAL)
                frq = np.arange(20, 201, 1)
                ice = 0
                gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
                rh = mr2rh(p, tk, gkg)[0] / 100

                e, rho = RTEquation.vapor(tk, rh, ice)

                AbsModel.model = 'R16'
                H2OAbsModel.h2oll = import_lineshape('h2oll')
                for i in range(0, len(z)):
                    v = 300.0 / tk[i]
                    ekpa = e[i] / 10.0
                    pdrykpa = p[i] / 10.0 - ekpa
                    for j in range(0, len(frq)):
                        _, _ = H2OAbsModel().h2o_absorption(pdrykpa, v, ekpa, frq[j])

        """
        if amu:
            if H2OAbsModel.model > 'R17':
                self.h2oll.cf = amu['con_Cf'].value
                self.h2oll.cs = amu['con_Cs'].value
            else:
                self.h2oll.cf = amu['con_Cf'].value * amu['con_Cf_factr'].value
                self.h2oll.cs = amu['con_Cs'].value * amu['con_Cs_factr'].value
            self.h2oll.xcf = amu['con_Xf'].value
            self.h2oll.xcs = amu['con_Xs'].value
            self.h2oll.s1 = amu['S'].value
            self.h2oll.b2 = amu['B2'].value
            self.h2oll.w0 = amu['gamma_a'].value / 1000.0
            self.h2oll.w0s = amu['gamma_w'].value / 1000.0
            self.h2oll.x = amu['n_a'].value
            self.h2oll.xs = amu['n_w'].value
            if H2OAbsModel.model > 'R17':
                self.h2oll.sh = amu['delta_a'].value / 1000.0
                self.h2oll.shs = amu['delta_w'].value / 1000.0
                self.h2oll.xh = amu['n_da'].value
                self.h2oll.xhs = amu['n_dw'].value
            else:
                self.h2oll.sr = amu['SR'].value
            self.h2oll.fl = amu['FL'].value

        # the best-fit voigt are given in koshelev et al. 2018, table 2 (rad,
        # mhz/torr). these correspond to w3(1) and ws(1) in h2o_list_r18 (mhz/mb)

        # cyh ***********************************************
        db2np = np.log(10.0) * 0.1
        rvap = (0.01 * 8.31451) / 18.01528
        factor = 0.182 * frq
        t = 300.0 / vx
        p = (pdrykpa + ekpa) * 10.0
        rho = ekpa * 10.0 / (rvap * t)
        f = frq
        # cyh ***********************************************

        if rho.any() <= 0.0:
            npp = 0
            ncpp = 0
            return npp, ncpp

        pvap = (rho * t) / 216.68
        if H2OAbsModel.model in ['R03', 'R16', 'R17', 'R98']:
            pvap = (rho * t) / 217.0
        if H2OAbsModel.model in ['R22SD', 'R23SD', 'R24']:
            pvap = (constants("Rwatvap")[0] * 1e-05) * rho * t
        pda = p - pvap
        if H2OAbsModel.model in ['R03', 'R16', 'R98']:
            den = 3.335e+16 * rho
        else:
            den = 3.344e+16 * rho
        # continuum terms
        ti = self.h2oll.reftcon / t
        # xcf and xcs include 3 for conv. to density & stimulated emission
        if H2OAbsModel.model in ['R03', 'R98']:
            con = (5.43e-10 * pda * ti ** 3 + 1.8e-08 *
                   pvap * ti ** 7.5) * pvap * f * f
        else:
            con = (self.h2oll.cf * pda * ti ** self.h2oll.xcf + self.h2oll.cs * pvap * ti ** self.h2oll.xcs) * \
                pvap * f * f
        # 2019/03/18 *********************************************************
        # add resonances
        nlines = len(self.h2oll.fl)
        ti = self.h2oll.reftline / t
        df = np.zeros((2, 1))

        if H2OAbsModel.model.startswith(('R19SD', 'R20SD', 'R21SD', 'R22SD', 'R23SD', 'R24')):
            tiln = np.log(ti)
            ti2 = np.exp(2.5 * tiln)
            summ = 0.0
            if H2OAbsModel.model in ['R23SD', 'R24']:
                if self.h2oll.cs > 0:
                    # npp_cs = np.zeros(1)
                    con = self.h2oll.cs * ti * self.h2oll.xcs
                    # for i in len(frq):
                    npp_cs = con
                else:
                    npp_cs = self.h2o_continuum(frq, vx, 1)
            for i in range(0, nlines):
                width0 = self.h2oll.w0[i] * pda * ti ** self.h2oll.x[i] + \
                    self.h2oll.w0s[i] * pvap * ti ** self.h2oll.xs[i]
                width2 = self.h2oll.w2[i] * pda + self.h2oll.w2s[i] * pvap
                if H2OAbsModel.model in ['R21SD', 'R22SD', 'R23SD', 'R24']:
                    if self.h2oll.w2[i] > 0:
                        width2 = self.h2oll.w2[i] * pda * ti ** self.h2oll.xw2[i] + self.h2oll.w2s[i] * pvap * ti ** \
                            self.h2oll.xw2s[i]
                    else:
                        width2 = 0

                shiftf = self.h2oll.sh[i] * pda * \
                    (1. - self.h2oll.aair[i] * tiln) * ti ** self.h2oll.xh[i]
                shifts = self.h2oll.shs[i] * pvap * \
                    (1. - self.h2oll.aself[i] * tiln) * ti ** self.h2oll.xhs[i]
                shift = shiftf + shifts
                # thus using the best-fit voigt (shift instead of shift0 and shift2)
                wsq = width0 ** 2
                s = self.h2oll.s1[i] * ti2 * \
                    np.exp(self.h2oll.b2[i] * (1. - ti))
                df[0] = f - self.h2oll.fl[i] - shift
                df[1] = f + self.h2oll.fl[i] + shift
                base = width0 / (562500.0 + wsq)
                if H2OAbsModel.model in ["R21SD", 'R22SD', 'R23SD', 'R24']:
                    delta2 = self.h2oll.d2[i] * pda + self.h2oll.d2s[i] * pvap
                res = 0.0
                for j in range(0, 2):
                    if width2 > 0 and j == 0 and np.abs(df[j]) < (10 * width0):
                        # speed-dependent resonant shape factor, minus base
                        xc = complex(
                            (width0 - np.dot(1.5, width2)), df[j]) / width2
                        if H2OAbsModel.model == 'R20SD':
                            if i == 1:
                                delta2 = (self.h2oll.d2air * pda) + \
                                    (self.h2oll.d2self * pvap)
                            else:
                                delta2 = 0.0
                            xc = complex((width0 - np.dot(1.5, width2)), df[j] + np.dot(1.5, delta2)) / complex(
                                width2, -delta2)
                        elif H2OAbsModel.model in ["R21SD", 'R22SD', 'R23SD', 'R24']:
                            xc = complex(
                                (width0 - 1.5 * width2), df[j] + 1.5 * delta2) / complex(width2, -delta2)

                        xrt = np.sqrt(xc)
                        pxw = 1.77245385090551603 * xrt * \
                            _dcerror(-np.imag(xrt), np.real(xrt))
                        sd = 2.0 * (1.0 - pxw) / (
                            width2 if H2OAbsModel.model not in ['R20SD', 'R21SD', 'R22SD', 'R23SD', 'R24'] else complex(width2, -delta2))
                        res += np.real(sd) - base
                    elif np.abs(df[j]) < 750.0:
                        res += width0 / (df[j] ** 2 + wsq) - base

                summ += s * res * (f / self.h2oll.fl[i]) ** 2
        elif H2OAbsModel.model in ['R16', 'R03', 'R17', 'R98']:
            ti2 = ti ** 2.5
            summ = 0.0
            for i in range(0, nlines):
                widthf = self.h2oll.w0[i] * pda * ti ** self.h2oll.x[i]
                widths = self.h2oll.w0s[i] * pvap * ti ** self.h2oll.xs[i]
                width = widthf + widths
                if H2OAbsModel.model == 'R98':
                    shift = 0.0
                else:
                    shift = self.h2oll.sr[i] * \
                        (width if H2OAbsModel.model == 'R03' else widthf)
                wsq = width ** 2
                s = self.h2oll.s1[i] * ti2 * \
                    np.exp(self.h2oll.b2[i] * (1.0 - ti))
                df[0] = f - self.h2oll.fl[i] - shift
                df[1] = f + self.h2oll.fl[i] + shift
                # use clough's definition of local line contribution
                base = width / (562500.0 + wsq)
                # do for positive and negative resonances
                res = 0.0
                for j in range(0, 2):
                    if np.abs(df[j]) <= 750.0:
                        res += width / (df[j] ** 2 + wsq) - base
                summ += s * res * (f / self.h2oll.fl[i]) ** 2
        elif H2OAbsModel.model in ['R19', 'R20']:
            tiln = np.log(ti)
            ti2 = ti ** 2.5
            summ = 0.0
            for i in range(0, nlines):
                widthf = self.h2oll.w0[i] * pda * ti ** self.h2oll.x[i]
                widths = self.h2oll.w0s[i] * pvap * ti ** self.h2oll.xs[i]
                width = widthf + widths
                shiftf = self.h2oll.sh[i] * pda * \
                    (1. - self.h2oll.aair[i] * tiln) * ti ** self.h2oll.xh[i]
                shifts = self.h2oll.shs[i] * pvap * \
                    (1. - self.h2oll.aself[i] * tiln) * ti ** self.h2oll.xhs[i]
                shift = shiftf + shifts
                wsq = width ** 2
                s = self.h2oll.s1[i] * ti2 * \
                    np.exp(self.h2oll.b2[i] * (1. - ti))
                df[0] = f - self.h2oll.fl[i] - shift
                df[1] = f + self.h2oll.fl[i] + shift
                base = width / (562500.0 + wsq)
                res = 0.0
                for j in range(0, 2):
                    if np.abs(df[j]) < 750.0:
                        res += width / (df[j] ** 2 + wsq) - base
                summ += s * res * (f / self.h2oll.fl[i]) ** 2
        elif H2OAbsModel.model == 'R18':
            ti2 = ti ** 2.5
            summ = 0.0
            for i in range(0, nlines):
                widthf = self.h2oll.w0[i] * pda * ti ** self.h2oll.x[i]
                widths = self.h2oll.w0s[i] * pvap * ti ** self.h2oll.xs[i]
                width = widthf + widths
                shiftf = self.h2oll.sh[i] * pda * ti ** self.h2oll.xh[i]
                shifts = self.h2oll.shs[i] * pvap * ti ** self.h2oll.xhs[i]
                shift = shiftf + shifts
                wsq = width ** 2
                s = self.h2oll.s1[i] * ti2 * \
                    np.exp(self.h2oll.b2[i] * (1. - ti))
                df[0] = f - self.h2oll.fl[i] - shift
                df[1] = f + self.h2oll.fl[i] + shift
                base = width / (562500.0 + wsq)
                res = 0.0
                for j in range(0, 2):
                    if np.abs(df[j]) < 750.0:
                        res += width / (df[j] ** 2 + wsq) - base
                summ += s * res * (f / self.h2oll.fl[i]) ** 2
        # separate the following original equ. into line and continuum
        # terms, and change the units from np/km to ppm
        # abh2o = .3183e-4*den*sum + con

        if H2OAbsModel.model in ['R23SD', 'R24']:
            conf = self.h2oll.cf * ti**self.h2oll.xcf
            con = (conf * pda + npp_cs * pvap) * pvap * f**2

        h20m = 2.9915075E-23  # mass of water molecule (g)
        if H2OAbsModel.model in ['R22SD', 'R23SD']:
            npp = 1.e-10 * rho * summ / (np.pi * h20m) / db2np / factor
        elif H2OAbsModel.model == 'R24':
            npp = 1.e-10 * (rho/h20m) * (summ/np.pi) / db2np / factor
        else:
            npp = (3.1831e-05 * den * summ / db2np) / factor

        ncpp = (con / db2np) / factor

        return npp, ncpp


class O2AbsModel(AbsModel):
    """This class contains the :math:`O_2` absorption model used in pyrtlib.
    """

    def __init__(self) -> None:
        super(O2AbsModel, self).__init__()
        self._o2ll = None

    @property
    def o2ll(self) -> types.MethodType:
        """Getter for line list coefficient"""
        return self._o2ll

    @o2ll.setter
    def o2ll(self, o2ll) -> None:
        """Setter for line list coefficient"""
        if o2ll and isinstance(o2ll, types.MethodType):
            self._o2ll = o2ll
        else:
            raise ValueError("Please enter a valid absorption model")

    @staticmethod
    def set_ll() -> None:
        if O2AbsModel.model not in O2AbsModel.implemented_models()['Oxygen']:
            raise AbsModelError(O2AbsModel.model,
                                f"Model {O2AbsModel.model} is not available. It is necessary to define oxygen absorption model manually")
        O2AbsModel.o2ll = import_lineshape("o2ll")

    def o2_absorption(self, pdrykpa: float, vx: float, ekpa: float, frq: float, amu: Optional[dict] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Returns power absorption coefficient due to oxygen in air in nepers/km.

        History:

        * 5/1/95  P. Rosenkranz
        * 11/5/97  P. Rosenkranz - 1- line modification.
        * 12/16/98 pwr - updated submm freq's and intensities from HITRAN96
        * 8/21/02  pwr - revised width at 425
        * 3/20/03  pwr - 1- line mixing and width revised
        * 9/29/04  pwr - new widths and mixing, using HITRAN intensities for all lines
        * 6/12/06  pwr - chg. T dependence of 1- line to 0.8
        * 10/14/08 pwr - moved isotope abundance back into intensities, added selected O16O18 lines.
        * 5/30/09  pwr - remove common block, add weak lines.
        * 12/18/14 pwr - adjust line broadening due to water vapor.
        * 9/29/18  pwr - 2nd-order line mixing
        * 8/20/19  pwr - adjust intensities according to Koshelev meas.

        Args:
            pdrykpa (numpy.ndarray): Dry air pressure (kPa).
            vx (numpy.ndarray): Theta (adim) - (normalised temperature 300/t(K)).
            ekpa (numpy.ndarray): Water vapor partial pressure (kPa).
            frq (numpy.ndarray): Frequency (GHz) - (valid 0-1000 GHz).

        Returns:
            [numpy.ndarray]: Oxigen line and continuum absorption terms (ppm)

        References
        ----------
        .. [1] P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING BY MICROWAVE RADIOMETRY 
            (M.A. Janssen, ed., 1993)(http://hdl.handle.net/1721.1/68611).
        .. [2] G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, pp.282-287 (2003).
        .. [3] M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
        .. [4] M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
        .. [5] B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
        .. [6] D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
        .. [7] D.S. Makarov et al, JQSRT v.112, pp.1420-1428 (2011).
        .. [8] M.A. Koshelev et al, JQSRT, v.154, pp.24-27 (2015).
        .. [9] M.A. Koshelev et al, JQSRT, v.169, pp.91-95 (2016).
        .. [10] M.A. Koshelev et al, JQSRT, v.196, pp.78?86 (2017).
        .. [11] D.S. Makarov et al, JQSRT v.243 (March 2020) doi:10.1016/j.jqsrt.2019.106798

        Line intensities from HITRAN2004.
        Non-resonant intensity from JPL catalog.

        .. note::
         1. The mm line-width coefficients are from Tretyakov et al 2005,
            Makarov et al 2008, and Koshelev et al 2016;
            submm line-widths from Golubiatnikov & Krupnov, except 234-GHz line width from Drouin.
            Mixing coeff. from Makarov's 2018 revision.
         2. The same temperature dependence (X) is used for submillimeter
            line widths as in the 60 GHz band: (1/T)**X (Koshelev et al 2016).
        """

        if amu:
            self.o2ll.w2a = amu['w2a'].value
            # self.o2ll.apu = amu['APU'].value
            self.o2ll.w300[0:38] = amu['O2gamma'].value[0:38]
            self.o2ll.w300[34:49] = amu['O2gamma_NL'].value[34:49]
            self.o2ll.wb300 = amu['WB300'].value
            self.o2ll.f = amu['O2FL'].value
            self.o2ll.s300 = amu['O2S'].value
            self.o2ll.be = amu['O2BE'].value
            self.o2ll.x11 = amu['X11'].value
            self.o2ll.x16 = amu['X16'].value
            self.o2ll.x05 = amu['X05'].value
            self.o2ll.x = amu['X05'].value
            self.o2ll.snr = amu['Snr'].value
            self.o2ll.ns = amu['O2_nS'].value
            if O2AbsModel.model > 'R19':
                self.o2ll.y0 = amu['y0'].value
                self.o2ll.y1 = amu['y1'].value
                self.o2ll.dnu0 = amu['dnu0'].value
                self.o2ll.dnu1 = amu['dnu1'].value
            else:
                self.o2ll.y300 = amu['Y300'].value
                self.o2ll.y300[34:49] = amu['Y300_NL'].value[34:49]
                self.o2ll.v = amu['O2_V'].value
                self.o2ll.v[34:49] = amu['O2_V_NL'].value[34:49]

        # *** add the following lines *************************
        db2np = np.log(10.0) * 0.1
        rvap = (0.01 * 8.314510) / 18.01528
        factor = 0.182 * frq
        temp = 300.0 / vx
        pres = (pdrykpa + ekpa) * 10.0
        vapden = (ekpa * 10.0) / (rvap * temp)
        freq = frq
        # *****************************************************

        th = 300.0 / temp
        th1 = th - 1.0
        b = th**self.o2ll.x
        preswv = vapden * temp / 216.68
        if O2AbsModel.model in ['R03', 'R16', 'R17', 'R18', 'R98']:
            preswv = vapden * temp / 217.0
        if O2AbsModel.model in ['R22', 'R22SD', 'R23', 'R23SD', 'R24']:
            preswv = 4.615228e-3 * vapden * temp
        presda = pres - preswv
        den = .001 * (presda * b + 1.2 * preswv * th)
        if O2AbsModel.model in ['R03', 'R16', 'R98']:
            den = 0.001 * (presda * b + 1.1 * preswv * th)
        if O2AbsModel.model == 'R03':
            den = 0.001 * (presda * th ** 0.9 + 1.1 * preswv * th)
        if O2AbsModel.model == 'R98':
            den = 0.001 * (presda + 1.1 * preswv) * th
        pe2 = den * den
        if O2AbsModel.model in ['R23', 'R23SD']:
            pe1 = .99 * den
            pe2 = pe1**2
        if O2AbsModel.model in ['R24']:
            pe1 = den  # [8] includes common TEMP-dependence
            pe2 = pe1**2
        dfnr = self.o2ll.wb300 * den

        # intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
        summ = 1.584e-17 * freq * freq * dfnr / \
            (th * (freq * freq + dfnr * dfnr))
        if O2AbsModel.model in ['R03', 'R16', 'R17', 'R18', 'R98']:
            summ = 0.0
        nlines = len(self.o2ll.f)
        if O2AbsModel.model in ['R03', 'R98']:
            for k in range(0, nlines):
                if k == 0:
                    df = self.o2ll.w300[0] * den
                else:
                    df = self.o2ll.w300[k] * den
                y = 0.001 * pres * b * \
                    (self.o2ll.y300[k] + self.o2ll.v[k] * th1)
                strr = self.o2ll.s300[k] * np.exp(-self.o2ll.be[k] * th1)
                sf1 = (df + (freq - self.o2ll.f[k]) * y) / \
                    ((freq - self.o2ll.f[k]) ** 2 + df * df)
                sf2 = (df - (freq + self.o2ll.f[k]) * y) / \
                    ((freq + self.o2ll.f[k]) ** 2 + df * df)
                summ += strr * (sf1 + sf2) * (freq / self.o2ll.f[k]) ** 2
        elif O2AbsModel.model in ['R17', 'R18', 'R19', 'R19SD']:
            for k in range(0, nlines):
                df = self.o2ll.w300[k] * den
                fcen = self.o2ll.f[k]
                y = den * (self.o2ll.y300[k] + self.o2ll.v[k] * th1)
                strr = self.o2ll.s300[k] * np.exp(-self.o2ll.be[k] * th1)
                sf1 = (df + (freq - fcen) * y) / \
                    ((freq - self.o2ll.f[k]) ** 2 + df * df)
                sf2 = (df - (freq + fcen) * y) / \
                    ((freq + self.o2ll.f[k]) ** 2 + df * df)
                summ += strr * (sf1 + sf2) * (freq / self.o2ll.f[k]) ** 2
        elif O2AbsModel.model in ['R23']:
            summ = 0.
            anorm = 0.
            a = np.zeros(nlines)
            g = np.zeros(nlines)
            for k in range(0, nlines):
                a[k] = self.o2ll.s300[k] * \
                    np.exp(-self.o2ll.be[k] * th1)/self.o2ll.f[k]**2
                g[k] = self.o2ll.g0[k] + self.o2ll.g1[k]*th1
                if k > 0 and k <= 37:
                    summ += a[k]*g[k]
                    anorm += a[k]
            off = summ/anorm
            summ = (1.584e-17/th) * dfnr / (freq * freq + dfnr * dfnr)
            for k in range(0, nlines):
                width = self.o2ll.w300[k] * den
                y = pe1*(self.o2ll.y300[k]+self.o2ll.y1[k]*th1)
                if k > 0 and k <= 37:
                    gfac = 1. + pe2 * (g[k] - off)
                else:
                    gfac = 1.
                fcen = self.o2ll.f[k] + pe2 * \
                    (self.o2ll.dnu0[k] + self.o2ll.dnu1[k] * th1)
                if k == 0 and np.abs(freq-fcen) < 10. * width:
                    width2 = .076 * width
                    xc = complex(width-1.5*width2, freq-fcen)/width2
                    xrt = np.sqrt(xc)
                    pxw = 1.77245385090551603 * xrt * \
                        _dcerror(-np.imag(xrt), np.real(xrt))
                    asd = complex(1., y)*2.*(1.-pxw)/width2
                    sf1 = np.real(asd)
                else:
                    sf1 = (width*gfac + (freq-fcen)*y) / \
                        ((freq-fcen)**2 + width**2)
                sf2 = (width*gfac - (freq+fcen)*y)/((freq+fcen)**2 + width**2)
                summ += a[k] * (sf1+sf2)
                if k == 37:
                    summ = np.maximum(summ, 0.)
        elif O2AbsModel.model in ['R24']:
            anorm = 0.
            wnr = self.o2ll.wb300 * pe1
            sumy = 1.584e-17 * self.o2ll.wb300
            sumg = 0.
            asq = 0.
            adjy = .99  # adjustment following update of line intensities
            a = np.zeros(nlines)
            y = np.zeros(nlines)
            g = np.zeros(nlines)
            for k in range(0, nlines):
                a[k] = self.o2ll.s300[k] * \
                    np.exp(-self.o2ll.be[k] * th1) * th/self.o2ll.f[k]**2
                y[k] = adjy * (self.o2ll.y300[k] + self.o2ll.y1[k]*th1)
                if k <= 37:
                    sumy += 2. * a[k] * \
                        (self.o2ll.w300[k] + y[k] * self.o2ll.f[k])
                g[k] = self.o2ll.g0[k] + self.o2ll.g1[k] * th1
                if k > 0 and k <= 37:
                    sumg += a[k] * g[k]
                    asq += a[k]**2
                    anorm += a[k]
            # The bias adjustment for Y is applied
            # from K=2 to 38; therefore it is normalized by ANORM summed from 2 to 38.
            sumy2 = sumy/(2. * anorm)
            ratio = sumg/asq
            for k in range(1, 38):
                y[k] -= sumy2/self.o2ll.f[k]  # bias adjustment
                g[k] -= a[k] * ratio  # makes G orthogonal to A

            summ = 1.584e-17 * wnr/(freq * freq + wnr * wnr)
            for k in range(0, nlines):
                width = self.o2ll.w300[k] * pe1
                yk = pe1 * y[k]
                if k > 0 and k <= 37:
                    gfac = 1. + pe2 * g[k]
                else:
                    gfac = 1.
                fcen = self.o2ll.f[k] + pe2 * \
                    (self.o2ll.dnu0[k] + self.o2ll.dnu1[k] * th1)
                if k == 0 and np.abs(freq-fcen) < 10. * width:
                    width2 = .076 * width
                    xc = complex(width - 1.5*width2, (freq-fcen))/width2
                    xrt = np.sqrt(xc)
                    pxw = 1.77245385090551603 * xrt * \
                        _dcerror(-np.imag(xrt), np.real(xrt))
                    asd = complex(1., yk) * 2. * (1.-pxw)/width2
                    sf1 = np.real(asd)
                else:
                    sf1 = (width*gfac + (freq-fcen)*yk) / \
                        ((freq-fcen)**2 + width**2)
                sf2 = (width*gfac - (freq+fcen)*yk)/((freq+fcen)**2 + width**2)
                summ += a[k] * (sf1+sf2)
        else:
            for k in range(0, nlines):
                df = self.o2ll.w300[k] * den
                y = den * (self.o2ll.y0[k] + self.o2ll.y1[k] * th1)
                dnu = pe2 * (self.o2ll.dnu0[k] + self.o2ll.dnu1[k] * th1)
                strr = self.o2ll.s300[k] * np.exp(-self.o2ll.be[k] * th1)
                del1 = freq - self.o2ll.f[k] - dnu
                del2 = freq + self.o2ll.f[k] + dnu
                gfac = 1. + pe2 * (self.o2ll.g0[k] + self.o2ll.g1[k] * th1)
                d1 = del1 * del1 + df * df
                d2 = del2 * del2 + df * df
                sf1 = (df * gfac + del1 * y) / d1
                sf2 = (df * gfac - del2 * y) / d2
                summ += strr * (sf1 + sf2) * (freq / self.o2ll.f[k]) ** 2

        o2abs = 1.6097e+11 * summ * presda * th ** 3
        if O2AbsModel.model in ['R23']:
            o2abs = 1.6097e+11 * summ * presda * freq * freq * th**3
        if O2AbsModel.model in ['R24']:
            o2abs = 1.6097e+11 * summ * presda * (freq * th)**2
        if O2AbsModel.model in ['R03', 'R98']:
            o2abs = 5.034e+11 * summ * presda * th ** 3 / 3.14159
        # o2abs = 1.004 * np.maximum(o2abs, 0.0)
        if O2AbsModel.model != 'R98':
            o2abs = np.maximum(o2abs, 0.0)
        if O2AbsModel.model in ['R20', 'R20SD', 'R22', 'R22SD']:
            o2abs = 1.004 * np.maximum(o2abs, 0.0)

        # *** ********************************************************
        # separate the equ. into line and continuum
        # terms, and change the units from np/km to ppm

        # intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17

        ncpp = 1.584e-17 * freq * freq * dfnr / \
            (th * (freq * freq + dfnr * dfnr))
        #  .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
        # a/(pi*k*t_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont
        # for units conversion (pa->hpa, hz->ghz)
        # pa2hpa=1e-2; hz2ghz=1e-9; m2cm=1e2; m2km=1e-3; pa2hpa^-1 * hz2ghz * m2cm^-2 * m2km^-1 = 1e-8
        # th^3 = th(from ideal gas law 2.13) * th(from the mw approx of stimulated emission 2.16 vs. 2.14) *
        # th(from the partition sum 2.20)
        if O2AbsModel.model in ['R03', 'R98']:
            ncpp = 1.6e-17 * freq * freq * dfnr / \
                (th * (freq * freq + dfnr * dfnr))
            ncpp *= 5.034e+11 * presda * th ** 3 / 3.14159
        else:
            ncpp *= 1.6097e11 * presda * th ** 3  # n/pi*sum0
        # change the units from np/km to ppm
        npp = (o2abs / db2np) / factor
        ncpp = (ncpp / db2np) / factor
        ncpp = 0 if O2AbsModel.model in [
            'R19', 'R19SD', 'R20', 'R20SD', 'R22', 'R22SD', 'R23', 'R24'] else ncpp

        return npp, ncpp


class O3AbsModel(AbsModel):
    """This class contains the :math:`O_3` absorption model used in pyrtlib.
    """

    def __init__(self) -> None:
        super(O3AbsModel, self).__init__()
        self._o3ll = None

    @property
    def o3ll(self) -> types.MethodType:
        """Getter for line list coefficient"""
        return self._o3ll

    @o3ll.setter
    def o3ll(self, o3ll) -> None:
        """Setter for line list coefficient"""
        if o3ll and isinstance(o3ll, types.MethodType):
            self._o3ll = o3ll
        else:
            raise ValueError("Please enter a valid absorption model")

    @staticmethod
    def set_ll() -> None:
        if O3AbsModel.model not in O3AbsModel.implemented_models()['Ozone']:
            raise AbsModelError(O3AbsModel.model,
                                f"Model {O3AbsModel.model} is not available. It is necessary to define ozone absorption model manually")
        O3AbsModel.o3ll = import_lineshape("o3ll")

    def o3_absorption(self, t: np.ndarray, p: np.ndarray, f: np.ndarray, o3n: np.ndarray, amu: Optional[dict] = None) -> np.ndarray:
        """This function computes power absorption coeff (Np/km) in the atmosphere 
        due to selcted lines of ozone (:math:`O_3`).

        Args:
            t (np.ndarray): Temperature (K)
            p (np.ndarray): Total pressure (mb)
            f (np.ndarray): Frequency (GHz)
            o3n (np.ndarray): Ozone number density (molecules/m^3)

        Returns:
            np.ndarray: Ozone power absorption coeff. (Np/km)
        """

        if amu:
            self.o3ll.fl = amu['O3_FL'].value
            self.o3ll.s1 = amu['O3_S1'].value
            self.o3ll.b = amu['O3_B'].value
            self.o3ll.w = amu['O3_W'].value
            self.o3ll.x = amu['O3_X'].value
            self.o3ll.sr = amu['O3_SR'].value

        if o3n.any() <= 0:
            abs_o3 = 0
            return abs_o3

        den = 1e-06 * o3n
        ti = self.o3ll.reftline / t
        ti2 = ti ** 2.5
        qvinv = 1.0 - np.exp(-1008.0 / t)
        # add resonances within 1 ghz of f.  most of the ozone is in the
        # stratosphere, so lines are relatively narrow, and lorentz shape
        # factor is ok.
        summ = 0.0
        nlines = len(self.o3ll.fl)
        if O3AbsModel.model in ["R22", "R22SD"]:
            for k in range(0, nlines):
                if self.o3ll.fl[k] > (f + 1.0):
                    break
                if self.o3ll.fl[k] >= (f - 1.0):
                    widthc = self.o3ll.w[k] * p * ti ** self.o3ll.x[k]
                    betad = .62065e-7 * self.o3ll.fl[k] * np.sqrt(t)
                    arg1 = (self.o3ll.fl[k]-f)/betad
                    arg2 = widthc/betad
                    s = self.o3ll.s1[k] * np.exp(self.o3ll.b[k] * (1.0 - ti))
                    summ += s * np.real(_dcerror(arg1, arg2))/betad

            abs_o3 = .56419e-4 * summ * qvinv * ti2 * den
        else:
            for k in range(0, nlines):
                if self.o3ll.fl[k] > (f + 1.0):
                    break
                if self.o3ll.fl[k] >= (f - 1.0):
                    widthc = self.o3ll.w[k] * p * ti ** self.o3ll.x[k]
                    betad2 = 3.85e-15 * t * self.o3ll.fl[k] ** 2
                    # approximate width combines pressure and doppler broadening:
                    width = 0.5346 * widthc + \
                        np.sqrt(0.2166 * widthc * widthc + 0.6931 * betad2)
                    s = self.o3ll.s1[k] * np.exp(self.o3ll.b[k] * (1.0 - ti))
                    shape = (f / self.o3ll.fl[k]) ** 2 * width / \
                        ((f - self.o3ll.fl[k]) ** 2 + width * width)
                    summ += s * shape

            abs_o3 = 3.183e-05 * summ * qvinv * ti2 * den

        return abs_o3

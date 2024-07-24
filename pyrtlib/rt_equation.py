# -*- coding: utf-8 -*-
"""
This class contains the main Radiative Transfer Equation functions.
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

import warnings
from typing import Tuple, Optional, Union

import numpy as np

from .absorption_model import O2AbsModel, H2OAbsModel, N2AbsModel, LiqAbsModel, O3AbsModel
from .utils import constants, tk2b_mod


class RTEquation:
    """This class contains the main Radiative Transfer Equation functions.
    """

    _from_sat = False
    _emissivity = 1.0

    @staticmethod
    def vapor(t: np.ndarray, rh: np.ndarray, ice: Optional[bool] = False) -> Tuple[np.ndarray, np.ndarray]:
        """Compute saturation vapor pressure (es,in mb) over water or ice at
        temperature t (kelvins), using the Goff-Gratch formulation (List,1963).

        Args:
            t (numpy.ndarray): Temperature (K).
            rh (numpy.ndarray): Relative Humidity (fraction).
            ice (Optional[bool], optional): Switch to calculate saturation vapor pressure over
                        water only (True) or water and ice, depending on T. Defaults to False.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: 
            * e (np.ndarray): Vapor pressure (mb)
            * rho (np.ndarray): Vapor density (:math:`g/m^3`)
        """

        rvap = constants('Rwatvap')[0]

        rvap = np.dot(rvap, 1e-05)

        # if ( (tk > 263.16) | (ice==0) )
        #    # for water...
        #    y = 373.16 ./ tk;
        #    es = -7.90298 * (y-1.) + 5.02808 * log10(y) -...
        #          1.3816e-7 * (10 .^ (11.344 * (1. - (1./ y))) - 1.) +...
        #          8.1328e-3 * (10 .^ (-3.49149 * (y - 1.)) - 1.) +...
        #          log10(1013.246);
        # else
        #    # for ice...
        #    y = 273.16 ./ tk;
        #    es = -9.09718 * (y - 1.) - 3.56654 * log10(y) +...
        #          0.876793 * (1.- (1. ./ y)) + log10(6.1071);
        # end

        # over water...
        y = 373.16 / t
        es = np.dot(-7.90298, (y - 1.0)) + np.dot(5.02808, np.log10(y)) - np.dot(1.3816e-07, (
            10 ** (np.dot(11.344, (1.0 - (1.0 / y)))) - 1.0)) + np.dot(0.0081328, (
                10 ** (np.dot(-3.49149, (y - 1.0))) - 1.0)) + np.log10(1013.246)
        if ice:
            # over ice if tk < 263.16
            indx = np.nonzero(t < 263.16)
            y = 273.16 / t(indx)
            es[indx] = np.dot(-9.09718, (y - 1.0)) - np.dot(
                3.56654, np.log10(y)) + np.dot(0.876793, (1.0 - (1.0 / y))) + np.log10(6.1071)

        es = 10.0 ** es
        # Compute vapor pressure and vapor density.
        # The vapor density conversion follows the ideal gas law:
        # vapor pressure = vapor density * rvapor * tk

        e = np.multiply(rh, es)
        rho = e / (np.dot(rvap, t))

        return e, rho

    @staticmethod
    def bright(hvk: np.ndarray, boft: np.ndarray) -> np.ndarray:
        r"""Function to compute temperature from the modified Planck
        radiance (Planck function without the constants :math:`\frac{2h\nu^3}{c^2}`.

        .. math:: B_{\nu}(\nu,T) = \frac{1}{ e^{\frac{h\nu}{k_{B}T}}-1}\implies T_b = \frac{h\nu}{k_{B}}\times\frac{1}{\ln(1+\frac{1}{B_{\nu}(\nu,T)})}

        Args:
            hvk (numpy.ndarray): (Planck constant (J*S)] * [frequency (Hz)) / [Boltzmann constant (J/K)].
            boft (numpy.ndarray): Modified Planck radiance - (equation (4) from [Schroeder-Westwater-1991]_).

        Returns:
            numpy.ndarray: Temperature (K)
        """

        if boft != 0:
            Tb = hvk / np.log(1.0 + (1.0 / boft))
        else:
            Tb = 0

        return Tb

    @staticmethod
    def refractivity(p: np.ndarray, t: np.ndarray, e: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Computes profiles of wet refractivity, dry refractivity, refractive index. 
        Refractivity equations were taken from [Thayer-1974]_.

        These equations were intended for frequencies under 20 GHz

        Args:
            p (numpy.ndarray): Pressure profile (mb). 
            t (numpy.ndarray): Temperature profile (K).
            e (numpy.ndarray): Vapor pressure profile (mb).

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:

            * dryn (numpy.ndarray): Dry refractivity profile
            * wetn (numpy.ndarray): Wet refractivity profile
            * refindx (numpy.ndarray): Refractive index profile
        """

        nl = len(p)
        wetn = np.zeros(p.shape)
        dryn = np.zeros(p.shape)
        refindx = np.zeros(p.shape)

        for i in range(0, nl):
            # Calculate dry air pressure (pa) and celsius temperature (tc).
            pa = p[i] - e[i]
            tc = t[i] - 273.16
            tk2 = np.dot(t[i], t[i])
            tc2 = np.dot(tc, tc)
            rza = 1.0 + \
                np.dot(
                    pa, (np.dot(5.79e-07, (1.0 + 0.52 / t[i])) - np.dot(0.00094611, tc) / tk2))
            rzw = 1.0 + np.dot(np.dot(1650.0, (e[i] / (np.dot(t[i], tk2)))),
                               (1.0 - np.dot(0.01317, tc) + np.dot(0.000175, tc2) + np.dot(1.44e-06,
                                                                                           (np.dot(tc2, tc)))))
            wetn[i] = np.dot((np.dot(64.79, (e[i] / t[i])) +
                             np.dot((377600.0), (e[i] / tk2))), rzw)
            dryn[i] = np.dot(np.dot(77.6036, (pa / t[i])), rza)
            refindx[i] = 1.0 + np.dot((dryn[i] + wetn[i]), 1e-06)

        return dryn, wetn, refindx

    @staticmethod
    def ray_tracing(z: np.ndarray, refindx: np.ndarray, angle: float, z0: float) -> Union[np.ndarray,  None]:
        """Ray-tracing algorithm of Dutton, Thayer, and Westwater, rewritten for 
        readability & attempted documentation. Based on the technique shown in Radio Meteorology 
        by Bean and Dutton (Fig. 3.20 and surrounding text) [Bean-Dutton]_.

        Args:
            z (numpy.ndarray): Height profile (km above observation height, z0).
            refindx (numpy.ndarray): Refractive index profile.
            angle (float): Elevation angle (degrees).
            z0 (float): Observation height (km msl).

        Returns:
            numpy.ndarray: Array containing slant path length profiles (km)

        .. note::
            The algorithm assumes that x decays exponentially over each layer.
        """

        deg2rad = np.pi / 180
        re = constants('EarthRadius')[0]
        ds = np.zeros(z.shape)

        nl = len(z)
        # Check for refractive index values that will blow up calculations.
        for i in range(0, nl):
            if refindx[i] < 1:
                warnings.warn('ray_tracing: Negative rafractive index')
                return

        # If angle is close to 90 degrees, make ds a height difference profile.
        if (angle >= 89 and angle <= 91) or (angle >= -91 and angle <= -89):
            ds[0] = 0.0
            for i in range(1, nl):
                ds[i] = z[i] - z[i - 1]
            return ds

        # The rest of the subroutine applies only to angle other than 90 degrees.
        # Convert angle degrees to radians.  Initialize constant values.
        theta0 = np.dot(angle, deg2rad)
        rs = re + z[0] + z0
        costh0 = np.cos(theta0)
        sina = np.sin(np.dot(theta0, 0.5))
        a0 = np.dot(2.0, (sina ** 2))
        # Initialize lower boundary values for 1st layer.
        ds[0] = 0.0
        phil = 0.0
        taul = 0.0
        rl = re + z[0] + z0
        tanthl = np.tan(theta0)
        # Construct the slant path length profile.
        for i in range(1, nl):
            r = re + z[i] + z0
            if refindx[i] == refindx[i - 1] or refindx[i] == 1. or refindx[i - 1] == 1.:
                refbar = np.dot((refindx[i] + refindx[i - 1]), 0.5)
            else:
                refbar = 1.0 + (refindx[i - 1] - refindx[i]) / \
                    (np.log((refindx[i - 1] - 1.0) / (refindx[i] - 1.0)))
            argdth = z[i] / rs - \
                (np.dot((refindx[0] - refindx[i]), costh0) / refindx[i])
            argth = np.dot(0.5, (a0 + argdth)) / r
            if argth <= 0:
                warnings.warn(
                    'ray_tracing: Ducting at {} degrees'.format(angle))
                return ds
            # Compute d-theta for this layer.
            sint = np.sqrt(np.dot(r, argth))
            theta = np.dot(2.0, np.arcsin(sint))
            if (theta - np.dot(2.0, theta0)) <= 0.0:
                dendth = np.dot(np.dot(2.0, (sint + sina)),
                                np.cos(np.dot((theta + theta0), 0.25)))
                sind4 = (np.dot(0.5, argdth) - np.dot(z[i], argth)) / dendth
                dtheta = np.dot(4.0, np.arcsin(sind4))
                theta = theta0 + dtheta
            else:
                dtheta = theta - theta0
            # Compute d-tau for this layer (eq.3.71) and add to integral, tau.
            tanth = np.tan(theta)
            cthbar = np.dot(((1.0 / tanth) + (1.0 / tanthl)), 0.5)
            dtau = np.dot(cthbar, (refindx[i - 1] - refindx[i])) / refbar
            tau = taul + dtau
            phi = dtheta + tau
            ds[i] = np.sqrt(
                (z[i] - z[i - 1]) ** 2 + np.dot(np.dot(np.dot(4.0, r), rl), ((np.sin(np.dot((phi - phil), 0.5))) ** 2)))
            if dtau != 0.0:
                dtaua = np.abs(tau - taul)
                ds[i] = np.dot(
                    ds[i], (dtaua / (np.dot(2.0, np.sin(np.dot(dtaua, 0.5))))))
            # Make upper boundary into lower boundary for next layer.
            phil = np.copy(phi)
            taul = np.copy(tau)
            rl = np.copy(r)
            tanthl = np.copy(tanth)

        return np.asarray(ds)

    @staticmethod
    def exponential_integration(zeroflg: bool, x: np.ndarray, ds: np.ndarray, ibeg: int, iend: int, factor: float) -> Tuple[np.ndarray, np.ndarray]:
        """EXPonential INTegration: Integrate the profile in array x over the layers defined in
        array ds, saving the integrals over each layer.

        Args:
            zeroflg (bool): Flag to handle zero values (0:layer=0, 1:layer=avg).
            x (numpy.ndarray): Profile array.
            ds (numpy.ndarray): Array of layer depths (km).
            ibeg (int): Lower integration limit (profile level number).
            iend (int): Upper integration limit (profile level number).
            factor (float): Factor by which result is multiplied (e.g., unit change).

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: 
            * xds (numpy.ndarray): Array containing integrals over each layer ds
            * sxds (numpy.ndarray): Integral of x*ds over levels ibeg to iend
        """

        sxds = 0.0
        xds = np.zeros(ds.shape)
        # TODO: check index
        for i in range(ibeg, iend):
            # Check for negative x value. If found, output message and return.
            if x[i - 1] < 0.0 or x[i] < 0.0:
                warnings.warn('Error encountered in exponential_integration')
                return sxds, xds
                # Find a layer value for x in cases where integration algorithm fails.
            elif np.abs(x[i] - x[i - 1]) < 1e-09:
                xlayer = x[i]
            elif x[i - 1] == 0.0 or x[i] == 0.0:
                if not zeroflg:
                    xlayer = 0.0
                else:
                    xlayer = np.dot((x[i] + x[i - 1]), 0.5)
            else:
                # Find a layer value for x assuming exponential decay over the layer.
                xlayer = (x[i] - x[i - 1]) / np.log(x[i] / x[i - 1])
            # Integrate x over the layer and save the result in xds.
            xds[i] = np.dot(xlayer, ds[i])
            sxds = sxds + xds[i]

        sxds = np.dot(sxds, factor)

        # TODO: reashape xds array
        return sxds, xds.reshape(iend)

    @staticmethod
    def cloud_radiating_temperature(ibase: float, itop: float, hvk: np.ndarray, tauprof: np.ndarray, boftatm: np.ndarray) -> Union[np.ndarray, None]:
        """Computes the mean radiating temperature of a cloud with base and top at
        profile levels ibase and itop, respectively. The algorithm assumes that
        the input cloud is the lowest (or only) cloud layer observed.
        If absorption is not too big, compute tmr of lowest cloud layer (base
        at level ibase, top at level itop). Otherwise, set error flag and return.

        Args:
            ibase (float): Profile level at base of lowest cloud.
            itop (float): Profile level at top of lowest cloud.
            hvk (np.ndarray): (Planck constant * frequency) / Boltzmann constant.
            tauprof (np.ndarray): Integral profile of absorption (np; i = integral (0,i)).
            boftatm (np.ndarray): Integral profile of atmospheric planck radiance.

        Returns:
            np.ndarray | None: tmr of lowest cloud layer (k)

        .. note::
            This algorithm is not designed for multiple cloud layers

        .. note::
            hvk, tauprof, and boftatm can be obtained from subroutine :py:meth:`planck`.
        """

        # maximum absolute value for exponential function argument
        expmax = 125.0
        ibase = int(ibase)
        itop = int(itop)
        # check if absorption too large to exponentiate
        if tauprof[ibase] > expmax:
            warnings.warn(
                'from cloud_radiating_temperature: absorption too large to exponentiate for tmr of lowest cloud layer')
            return

        # compute radiance (batmcld) and absorption (taucld) for cloud layer.
        # (if taucld is too large to exponentiate, treat it as infinity.)
        batmcld = boftatm[itop] - boftatm[ibase]
        taucld = tauprof[itop] - tauprof[ibase]
        if taucld > expmax:
            boftcld = np.dot(batmcld, np.exp(tauprof[ibase]))
        else:
            boftcld = (
                batmcld * np.exp(tauprof[ibase])) / (1.0 - np.exp(-taucld))

        # compute cloud mean radiating temperature (tmrcld)
        tmrcld = RTEquation.bright(hvk, boftcld)

        return tmrcld

    @staticmethod
    def cloud_integrated_density(dencld: np.ndarray, ds: np.ndarray, lbase: np.ndarray, ltop: np.ndarray) -> np.ndarray:
        """Integrates cloud water density over path ds (linear algorithm).

        Args:
            dencld (np.ndarray): Cloud cloud water density profile (:math:`g/m^3`)
            ds (np.ndarray): Vector containing layer depth profiles (km)
            lbase (np.ndarray): Array containing profile levels corresponding to cloud bases.
            ltop (np.ndarray): Array containing profile levels corresponding to cloud tops.

        Returns:
            np.ndarray: integrated cloud water density (cm)
        """

        ncld = len(lbase)
        scld = 0.0
        for i in range(0, ncld):
            for j in range(int(lbase[i]) + 1, int(ltop[i])):
                scld += np.dot(ds[j],
                               (np.dot(0.5, (dencld[j] + dencld[j - 1]))))

        # convert the integrated value to cm.
        scld = np.dot(scld, 0.1)

        return scld

    @staticmethod
    def planck(frq: np.ndarray, t: np.ndarray, taulay: np.ndarray) -> Tuple[
            np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Computes the modified planck function (equation (4) in [Schroeder-Westwater-1992]_ 
        for the cosmic background temperature, the mean radiating
        temperature, and a profile of the atmospheric integral with and without
        the cosmic background. Also computes an integral profile of atmospheric
        absorption. For the integral profiles, the value at profile level i
        represents the integral from the antenna to level i.
        Also returns the cosmic background term for the rte.

        Args:
            frq (numpy.ndarray): Channel frequency (GHz).
            t (numpy.ndarray): Temperature profile (K).
            taulay (numpy.ndarray): Number of profile levels.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

            * hvk: (Planck constant * frequency) / Boltzmann constant.
            * boft: Modified planck function for raob temperature profile.
            * bakgrnd: Background term of radiative transfer equation.
            * boftatm: Array of atmospheric planck radiance integrated (0,i).
            * boftotl: Total planck radiance from the atmosphere plus bakgrnd.
            * boftmr: Modified planck function for mean radiating temperature.
            * tauprof: Array of integrated absorption (np; 0,i).
        """

        Tc = constants('Tcosmicbkg')[0]
        h = constants('planck')[0]
        k = constants('boltzmann')[0]
        fHz = np.dot(frq, 1e9)

        hvk = np.dot(fHz, h) / k
        # maximum absolute value for exponential function argument
        expmax = 125.0
        nl = len(t)
        tauprof = np.zeros(taulay.shape)
        boftatm = np.zeros(taulay.shape)
        boft = np.zeros(taulay.shape)

        if RTEquation._from_sat:
            boftotl = 0.0
            ###########################################################################
            # Then compute upwelling radiance
            # Adapted from Planck_xxx.m, but from Satellite i-1 becomes i+1
            # taulay changed to i+1, debugged by ISMAR project
            ###########################################################################
            Ts = t[0]
            boft[nl - 1] = tk2b_mod(hvk, t[nl - 1])
            for i in range(nl - 2, -1, -1):
                boft[i] = tk2b_mod(hvk, t[i])
                boftlay = (boft[i + 1] + np.dot(boft[i],
                           np.exp(-taulay[i+1]))) / (1.0 + np.exp(-taulay[i+1]))
                batmlay = np.dot(
                    np.dot(boftlay, np.exp(-tauprof[i + 1])), (1.0 - np.exp(-taulay[i+1])))
                boftatm[i] = boftatm[i + 1] + batmlay
                tauprof[i] = tauprof[i + 1] + taulay[i+1]

            # The background is a combination of surface emission and downwelling
            # radiance (boftotl) reflected by the surface
            if tauprof[0] < expmax:
                boftbg = np.dot(RTEquation._emissivity, tk2b_mod(
                    hvk, Ts)) + np.dot((1 - RTEquation._emissivity), boftotl)
                # boftbg_sat  = es * TK2B_mod(hvk,Ts); # SAT: eps * B(Tsrf) + (1-eps) B_dw
                bakgrnd = np.dot(boftbg, np.exp(-tauprof[0]))
                boftotl = bakgrnd + boftatm[0]
                boftmr = boftatm[0] / (1.0 - np.exp(-tauprof[0]))
            else:
                bakgrnd = 0.0
                boftotl = boftatm[0]
                boftmr = boftatm[0]
        else:
            boft[0] = tk2b_mod(hvk, t[0])
            for i in range(1, nl):
                boft[i] = tk2b_mod(hvk, t[i])
                boftlay = (boft[i - 1] + boft[i] *
                           np.exp(-taulay[i])) / (1.0 + np.exp(-taulay[i]))
                batmlay = boftlay * \
                    np.exp(- tauprof[i - 1]) * (1.0 - np.exp(-taulay[i]))
                boftatm[i] = boftatm[i - 1] + batmlay
                tauprof[i] = tauprof[i - 1] + taulay[i]
            # compute the cosmic background term of the rte; compute total planck
            # radiance for atmosphere and cosmic background; if absorption too large
            # to exponentiate, assume cosmic background was completely attenuated.
            if tauprof[nl - 1] < expmax:
                boftbg = tk2b_mod(hvk, Tc)
                bakgrnd = boftbg * np.exp(-tauprof[nl - 1])
                boftotl = bakgrnd + boftatm[nl - 1]
                boftmr = boftatm[nl - 1] / (1.0 - np.exp(-tauprof[nl - 1]))
            else:
                bakgrnd = 0.0
                boftotl = boftatm[nl - 1]
                boftmr = boftatm[nl - 1]

        return boftotl, boftatm, boftmr, tauprof, hvk, boft, bakgrnd

    @staticmethod
    def cloudy_absorption(t: np.ndarray, denl: np.ndarray, deni: np.ndarray, frq: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Multiplies cloud density profiles by a given fraction and computes the
        corresponding cloud liquid and ice absorption profiles, using Rosenkranz's
        cloud liquid absorption and ice absorption by [Westwater-1972]_.

        Args:
            t (numpy.ndarray): Temperature profile (k).
            denl (numpy.ndarray): Liquid density profile (:math:`g/m^3`).
            deni (numpy.ndarray): Ice density profile (:math:`g/m^3`).
            frq (numpy.ndarray): Frequency array (GHz).

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]: 
            * aliq: Liquid absorption profile (np/km)
            * aice: Ice absorption profile (np/km)

        See also:
            :py:func:`~pyrtlib.absorption_model.LiqAbsModel.liquid_water_absorption`

        """

        nl = len(t)
        c = np.dot(constants('light')[0], 100)

        ghz2hz = 1e9
        db2np = np.dot(np.log(10.0), 0.1)

        wave = c / (np.dot(frq, ghz2hz))

        aliq = np.zeros(denl.shape)
        aice = np.zeros(denl.shape)
        for i in range(0, nl):
            # Compute liquid absorption np/km.
            if denl[i] > 0:
                aliq[i] = LiqAbsModel.liquid_water_absorption(
                    denl[i], frq, t[i])
            # compute ice absorption (db/km); convert non-zero value to np/km.
            if deni[i] > 0:
                aice[i] = np.dot(
                    np.dot((8.18645 / wave), deni[i]), 0.000959553)
                aice[i] = np.dot(aice[i], db2np)

        return aliq, aice

    @staticmethod
    def clearsky_absorption(p: np.ndarray, t: np.ndarray, e: np.ndarray, frq: np.ndarray, o3n: Optional[np.ndarray] = None, amu: Optional[dict] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Computes profiles of water vapor and dry air absorption 
        for a given set of frequencies. Subroutines :math:`H_2O` and :math:`O_2` 
        contain the absorption model of [Liebe-Layton]_ with oxygen 
        interference coefficients from [Rosenkranz-1988]_.

        Args:
            p (numpy.ndarray): Pressure profile (mb).
            t (numpy.ndarray): Temperature profile (K).
            e (numpy.ndarray): Vapor pressure profile (mb).
            frq (numpy.ndarray): Frequency (GHz).
            o3n (numpy.ndarray, optional): Ozone Number Density (molecules/m3). Defaults to None.

        Raises:
            ValueError: Raises error if absorption model is not defined.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray]:
            * awet: Water vapor absorption profile (np/km).
            * adry: Dry air absorption profile (np/km).

        See also:
            :py:func:`~pyrtlib.absorption_model.H2OAbsModel.h2o_absorption`
            :py:func:`~pyrtlib.absorption_model.O2AbsModel.o2_absorption`
        """

        nl = len(p)
        awet = np.zeros(p.shape)
        adry = np.zeros(p.shape)
        aO2 = np.zeros(p.shape)
        aN2 = np.zeros(p.shape)
        aO3 = np.zeros(p.shape)
        factor = np.dot(0.182, frq)
        db2np = np.dot(np.log(10.0), 0.1)
        for i in range(0, nl):
            # Compute inverse temperature parameter; convert wet and dry p to kpa.
            v = 300.0 / t[i]
            ekpa = e[i] / 10.0
            pdrykpa = p[i] / 10.0 - ekpa
            # add H2O term
            npp, ncpp = H2OAbsModel().h2o_absorption(pdrykpa, v, ekpa, frq, amu)
            awet[i] = factor * (npp + ncpp) * db2np
            # add O2 term
            npp, ncpp = O2AbsModel().o2_absorption(pdrykpa, v, ekpa, frq, amu)
            aO2[i] = factor * (npp + ncpp) * db2np
            # add N2 term
            if N2AbsModel.model not in ['R03', 'R16', 'R17', 'R18', 'R98']:
                aN2[i] = N2AbsModel.n2_absorption(
                    t[i], np.dot(pdrykpa, 10), frq)

            if not N2AbsModel.model:
                raise ValueError(
                    'No model avalaible with this name: {} . Sorry...'.format('model'))

            if isinstance(o3n, np.ndarray) and O3AbsModel.model in ['R18', 'R21', 'R21SD', 'R22', 'R22SD', 'R23', 'R23SD', 'R24']:
                aO3[i] = O3AbsModel().o3_absorption(
                    t[i], p[i], frq, o3n[i], amu)

            adry[i] = aO2[i] + aN2[i] + aO3[i]

        return awet, adry

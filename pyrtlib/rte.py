# -*- coding: utf-8 -*-
"""
This class contains the main Radiative Transfer Equation functions.
"""
import warnings

import numpy as np

from .absmodel import O2AbsModel, H2OAbsModel, N2AbsModel, LiqAbsModel
from .utils import constants, tk2b_mod


class RTEquation:
    """This class contains the main Radiative Transfer Equation functions.
    """

    from_sat = False

    @staticmethod
    def vapor(tk=None, rh=None, ice=None, *args, **kwargs):
        """Compute saturation vapor pressure (es,in mb) over water or ice at
        temperature tk (kelvins), using the Goff-Gratch formulation (List,1963).

        Args:
            tk ([type], optional): temperature (K). Defaults to None.
            rh ([type], optional): relative humidity (fraction). Defaults to None.
            ice ([type], optional): switch to calculate saturation vapor pressure over
            water only (0) or water and ice, depending on tk (1). Defaults to None.

        Returns:
            (tuple):

            * e [type]: vapor pressure (mb)
            * rho [type]: vapor density (g/m3)

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
        y = 373.16 / tk
        es = np.dot(-7.90298, (y - 1.0)) + np.dot(5.02808, np.log10(y)) - np.dot(1.3816e-07, (
                10 ** (np.dot(11.344, (1.0 - (1.0 / y)))) - 1.0)) + np.dot(0.0081328, (
                10 ** (np.dot(-3.49149, (y - 1.0))) - 1.0)) + np.log10(1013.246)
        if ice == 1:
            # over ice if tk < 263.16
            indx = np.nonzero(tk < 263.16)
            y = 273.16 / tk(indx)
            es[indx] = np.dot(-9.09718, (y - 1.0)) - np.dot(
                3.56654, np.log10(y)) + np.dot(0.876793, (1.0 - (1.0 / y))) + np.log10(6.1071)

        es = 10.0 ** es
        # Compute vapor pressure and vapor density.
        # The vapor density conversion follows the ideal gas law:
        # apor pressure = vapor density * rvapor * tk

        e = np.multiply(rh, es)
        rho = e / (np.dot(rvap, tk))

        return e, rho

    @staticmethod
    def bright(hvk=None, boft=None, *args, **kwargs):
        """Function to compute temperature from the modified Planck
        radiance (Planck function without the constants 2h(v^3)/(c^2).

        Args:
            hvk ([type], optional): [Planck constant (J*S)] * [frequency (Hz)] / [Boltzmann constant (J/K)]. Defaults to None.
            boft ([type], optional): modified Planck radiance - (equation (4) from Schroeder & Westwater, 1991). Defaults to None.

        Returns:
            [type]: [description]
        """

        # fixes floa FloatingPointError: divide by zero encountered in double_scalars
        if boft != 0:
            Tb = hvk / np.log(1.0 + (1.0 / boft))
        else:
            Tb = np.nan

        return Tb

    @staticmethod
    def refractivity(p=None, tk=None, e=None, *args, **kwargs):
        """Computes profiles of wet refractivity, dry refractivity,
        refractive index.  Refractivity equations were taken from G.D.
        Thayer, 1974:  An improved equation for the radio refractive
        index of air. Radio Science, vol.9,no.10, 803-807.
        These equations were intended for frequencies under 20 GHz

        Args:
            p ([type], optional): pressure profile (mb). Defaults to None.
            tk ([type], optional): temperature profile (K). Defaults to None.
            e ([type], optional): vapor pressure profile (mb). Defaults to None.

        Returns:
            (tuple):

            * dryn [type]: dry refractivity profile
            * wetn [type]: wet refractivity profile
            * refindx [type]: refractive index profile

        .. todo:: check if slice works properly @Donatello
        """

        nl = len(p)
        wetn = np.zeros(p.shape)
        dryn = np.zeros(p.shape)
        refindx = np.zeros(p.shape)

        for i in range(0, nl):
            # Calculate dry air pressure (pa) and celsius temperature (tc).
            pa = p[i] - e[i]
            tc = tk[i] - 273.16
            tk2 = np.dot(tk[i], tk[i])
            tc2 = np.dot(tc, tc)
            rza = 1.0 + np.dot(pa, (np.dot(5.79e-07, (1.0 + 0.52 / tk[i])) - np.dot(0.00094611, tc) / tk2))
            rzw = 1.0 + np.dot(np.dot(1650.0, (e[i] / (np.dot(tk[i], tk2)))),
                               (1.0 - np.dot(0.01317, tc) + np.dot(0.000175, tc2) + np.dot(1.44e-06,
                                                                                           (np.dot(tc2, tc)))))
            wetn[i] = np.dot((np.dot(64.79, (e[i] / tk[i])) + np.dot((377600.0), (e[i] / tk2))), rzw)
            dryn[i] = np.dot(np.dot(77.6036, (pa / tk[i])), rza)
            refindx[i] = 1.0 + np.dot((dryn[i] + wetn[i]), 1e-06)

        return dryn, wetn, refindx

    @staticmethod
    def ray_tracing(z=None, refindx=None, angle=None, z0=None, *args, **kwargs):
        """Ray-tracing algorithm of Dutton, Thayer, and Westwater, rewritten for
        readability & attempted documentation.  Based on the technique shown in
        Radio Meteorology by Bean and Dutton (Fig. 3.20 and surrounding text).

        Args:
            z ([type], optional): eight profile (km above observation height, z0). Defaults to None.
            refindx ([type], optional): refractive index profile. Defaults to None.
            angle ([type], optional): elevation angle (degrees). Defaults to None.
            z0 ([type], optional): observation height (km msl). Defaults to None.

        Returns:
            [type]: array containing slant path length profiles (km)

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
                warnings.warn('RayTrac_xxx: Negative rafractive index')
                return

        # If angle is close to 90 degrees, make ds a height difference profile.
        if (angle >= 89 and angle <= 91) or (angle >= -91 and angle <= -89):
            ds[0] = 0.0
            for i in range(1, nl):
                ds[i] = z[i] - z[i - 1]

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
                refbar = 1.0 + (refindx[i - 1] - refindx[i]) / (np.log((refindx[i - 1] - 1.0) / (refindx[i] - 1.0)))
            argdth = z[i] / rs - (np.dot((refindx[0] - refindx[i]), costh0) / refindx[i])
            argth = np.dot(0.5, (a0 + argdth)) / r
            if argth <= 0:
                warnings.warn('RayTrac_xxx: Ducting at {} degrees'.format(angle))
                return ds
            # Compute d-theta for this layer.
            sint = np.sqrt(np.dot(r, argth))
            theta = np.dot(2.0, np.arcsin(sint))
            if (theta - np.dot(2.0, theta0)) <= 0.0:
                dendth = np.dot(np.dot(2.0, (sint + sina)), np.cos(np.dot((theta + theta0), 0.25)))
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
                ds[i] = np.dot(ds[i], (dtaua / (np.dot(2.0, np.sin(np.dot(dtaua, 0.5))))))
            # Make upper boundary into lower boundary for next layer.
            phil = np.copy(phi)
            taul = np.copy(tau)
            rl = np.copy(r)
            tanthl = np.copy(tanth)

            return np.asarray(ds)

    @staticmethod
    def exponential_integration(zeroflg=None, x=None, ds=None, ibeg=None, iend=None, factor=None, *args, **kwargs):
        """ EXPonential INTegration: Integrate the profile in array x over the layers defined in
        array ds, saving the integrals over each layer.

        Args:
            zeroflg ([type], optional): flag to handle zero values (0:layer=0, 1:layer=avg). Defaults to None.
            x ([type], optional): profile array. Defaults to None.
            ds ([type], optional): array of layer depths (km). Defaults to None.
            ibeg ([type], optional): lower integration limit (profile level number). Defaults to None.
            iend ([type], optional): upper integration limit (profile level number). Defaults to None.
            factor ([type], optional): factor by which result is multiplied (e.g., unit change). Defaults to None.

        Returns:
                (tuple):

                * xds [type]: array containing integrals over each layer ds
                * sxds [type]: integral of x*ds over levels ibeg to iend
        """

        sxds = 0.0
        xds = np.zeros(ds.shape)
        # TODO: check index
        for i in range(ibeg, iend):
            # Check for negative x value. If found, output message and return.
            if x[i - 1] < 0.0 or x[i] < 0.0:
                warnings.warn('Error encountered in ExpInt_xxx.m')
                return sxds, xds
                # Find a layer value for x in cases where integration algorithm fails.
            elif np.abs(x[i] - x[i - 1]) < 1e-09:
                xlayer = x[i]
            elif x[i - 1] == 0.0 or x[i] == 0.0:
                if zeroflg == 0:
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
    def cloud_radiating_temperature(ibase=None, itop=None, hvk=None, tauprof=None, boftatm=None, *args, **kwargs):
        """Computes the mean radiating temperature of a cloud with base and top at
        profile levels ibase and itop, respectively.  The algorithm assumes that
        the input cloud is the lowest (or only) cloud layer observed.
        If absorption is not too big, compute tmr of lowest cloud layer (base
        at level ibase, top at level itop). Otherwise, set error flag and return.

        Args:
            ibase ([type], optional): profile level at base of lowest cloud. Defaults to None.
            itop ([type], optional): profile level at top of lowest cloud. Defaults to None.
            hvk ([type], optional): (planck constant * frequency) / boltzmann constant. Defaults to None.
            tauprof ([type], optional): integral profile of absorption (np; i = integral (0,i)). Defaults to None.
            boftatm ([type], optional): integral profile of atmospheric planck radiance. Defaults to None.

        Returns:
            [type]: tmr of lowest cloud layer (k)

        .. note::
            This algorithm is not designed for multiple cloud layers

        .. note::
            hvk, tauprof, and boftatm can be obtained from subroutine planck_xxx().
        """

        # maximum absolute value for exponential function argument
        expmax = 125.0
        # ... check if absorption too large to exponentiate...
        if tauprof[ibase] > expmax:
            warnings.warn('from CldTmr_xxx: absorption too large to exponentiate for tmr of lowest cloud layer')
            return

        # compute radiance (batmcld) and absorption (taucld) for cloud layer.
        # (if taucld is too large to exponentiate, treat it as infinity.)
        batmcld = boftatm[itop] - boftatm[ibase]
        taucld = tauprof[itop] - tauprof[ibase]
        if taucld > expmax:
            boftcld = np.dot(batmcld, np.exp(tauprof[ibase]))
        else:
            boftcld = np.dot(batmcld, np.exp(tauprof[ibase])) / (1.0 - np.exp(-taucld))

        # compute cloud mean radiating temperature (tmrcld)
        tmrcld = RTEquation.bright(hvk, boftcld)

        return tmrcld

    @staticmethod
    def cloud_integrated_density(dencld=None, ds=None, lbase=None, ltop=None, *args, **kwargs):
        """Integrates cloud water density over path ds (linear algorithm).

        Args:
            dencld ([type], optional): cloud cloud water density profile (g/m3). Defaults to None.
            ds ([type], optional): vector containing layer depth profiles (km). Defaults to None.
            nlay ([type], optional): number of cloud layers in the profile. Defaults to None.
            lbase ([type], optional): array containing profile levels corresponding to cloud bases. Defaults to None.
            ltop ([type], optional): array containing profile levels corresponding to cloud tops . Defaults to None.

        Returns:
            [type]: integrated cloud water density (cm)

        .. warning:: nlay arg is not defined ??
        """
        ncld = len(lbase)
        scld = 0.0
        for i in range(0, ncld):
            for j in range(lbase[i] + 1, ltop[i]):
                scld = scld + np.dot(ds[j], (np.dot(0.5, (dencld[j] + dencld[j - 1]))))

        # convert the integrated value to cm.
        scld = np.dot(scld, 0.1)

        return scld

    @staticmethod
    def planck(frq=None, tk=None, taulay=None, *args, **kwargs):
        """  Computes the modified planck function (equation (4) in schroeder and
        westwater, 1992: guide to passive microwave weighting function
        calculations) for the cosmic background temperature, the mean radiating
        temperature, and a profile of the atmospheric integral with and without
        the cosmic background. Also computes an integral profile of atmospheric
        absorption. For the integral profiles, the value at profile level i
        represents the integral from the antenna to level i.
        Also returns the cosmic background term for the rte.

        Args:
            frq ([type], optional): channel frequency (GHz). Defaults to None.
            nl ([type], optional): number of profile levels
            tk ([type], optional): temperature profile (K). Defaults to None.
            taulay ([type], optional): profile of absorption integrated over each layer (np). Defaults to None.

        Returns:
            (tuple):

                * hvk [type]: [planck constant * frequency] / boltzmann constant
                * boft [type]: modified planck function for raob temperature profile
                * bakgrnd [type]: background term of radiative transfer equation
                * boftatm [type]: array of atmospheric planck radiance integrated (0,i)
                * boftotl [type]: total planck radiance from the atmosphere plus bakgrnd
                * boftmr  [type]: modified planck function for mean radiating temperature
                * tauprof [type]: array of integrated absorption (np; 0,i)

        .. warning:: nl arg is missing ??
        """

        Tc = constants('Tcosmicbkg')[0]
        h = constants('planck')[0]
        k = constants('boltzmann')[0]
        fHz = np.dot(frq, 1000000000.0)

        hvk = np.dot(fHz, h) / k
        # maximum absolute value for exponential function argument
        expmax = 125.0
        nl = len(tk)
        tauprof = np.zeros(taulay.shape)
        boftatm = np.zeros(taulay.shape)
        boft = np.zeros(taulay.shape)

        if RTEquation.from_sat:
            boftotl = 0.0
            ###########################################################################
            # Then compute upwelling radiance 
            # Adapted from Planck_xxx.m, but from Satellite i-1 becomes i+1
            ###########################################################################
            Ts = tk[0]
            Es = 1.0
            boft[nl - 1] = tk2b_mod(hvk, tk[nl - 1])
            for i in range(nl - 2, -1, -1):
                boft[i] = tk2b_mod(hvk, tk[i])
                boftlay = (boft[i + 1] + np.dot(boft[i], np.exp(-taulay[i]))) / (1.0 + np.exp(-taulay[i]))
                batmlay = np.dot(np.dot(boftlay, np.exp(-tauprof[i + 1])), (1.0 - np.exp(-taulay[i])))
                boftatm[i] = boftatm[i + 1] + batmlay
                tauprof[i] = tauprof[i + 1] + taulay[i]

            # The background is a combination of surface emission and downwelling 
            # radiance (boftotl) reflected by the surface
            if tauprof[0] < expmax:
                boftbg = np.dot(Es, tk2b_mod(hvk, Ts)) + np.dot((1 - Es), boftotl)
                # boftbg_sat  = Es * TK2B_mod(hvk,Ts); # SAT: eps * B(Tsrf) + (1-eps) B_dw
                bakgrnd = np.dot(boftbg, np.exp(-tauprof[0]))
                boftotl = bakgrnd + boftatm[0]
                boftmr = boftatm[0] / (1.0 - np.exp(-tauprof[0]))
            else:
                bakgrnd = 0.0
                boftotl = boftatm[0]
                boftmr = boftatm[0]
        else:
            # TODO: check index of boft variable
            boft[0] = tk2b_mod(hvk, tk[0])
            for i in range(1, nl):
                boft[i] = tk2b_mod(hvk, tk[i])
                boftlay = (boft[i - 1] + np.dot(boft[i], np.exp(-taulay[i]))) / (1.0 + np.exp(-taulay[i]))
                batmlay = np.dot(np.dot(boftlay, np.exp(- tauprof[i - 1])), (1.0 - np.exp(-taulay[i])))
                boftatm[i] = boftatm[i - 1] + batmlay
                tauprof[i] = tauprof[i - 1] + taulay[i]
            # compute the cosmic background term of the rte; compute total planck
            # radiance for atmosphere and cosmic background; if absorption too large
            # to exponentiate, assume cosmic background was completely attenuated.
            if tauprof[nl - 1] < expmax:
                boftbg = tk2b_mod(hvk, Tc)
                bakgrnd = np.dot(boftbg, np.exp(-tauprof[nl - 1]))
                boftotl = bakgrnd + boftatm[nl - 1]
                boftmr = boftatm[nl - 1] / (1.0 - np.exp(-tauprof[nl - 1]))
            else:
                bakgrnd = 0.0
                boftotl = boftatm[nl - 1]
                boftmr = boftatm[nl - 1]

        return boftotl, boftatm, boftmr, tauprof, hvk, boft, bakgrnd

    @staticmethod
    def cloudy_absorption(tk=None, denl=None, deni=None, frq=None, *args, **kwargs):
        """Multiplies cloud density profiles by a given fraction and computes the
        corresponding cloud liquid and ice absorption profiles, using Rosenkranz's
        cloud liquid absorption routine ABLIQ and ice absorption of Westwater
        [1972: Microwave Emission from Clouds,13-14]. - Yong Han, 4/20/2000

        Args:
            tk ([type], optional): temperature profile (k). Defaults to None.
            denl ([type], optional): liquid density profile (g/m3) . Defaults to None.
            deni ([type], optional): ice density profile (g/m3). Defaults to None.
            frq ([type], optional): frequency array (GHz). Defaults to None.

        Returns:
            (tuple):

                * aliq [type]: liquid absorption profile (np/km)
                * aice [type]: ice absorption profile (np/km)

        See also:

            :py:meth:`absmodel.AbsModel.ab_liq`

        .. warning::
            * ic light speed in cm s-1????
            * ab_liq function is missing!!!!!
        """

        nl = len(tk)
        c = np.dot(constants('light')[0], 100)

        ghz2hz = 1000000000.0
        db2np = np.dot(np.log(10.0), 0.1)

        wave = c / (np.dot(frq, ghz2hz))

        aliq = np.zeros(denl.shape)
        aice = np.zeros(denl.shape)
        for i in range(0, nl):
            # Compute liquid absorption np/km.
            if denl[i] > 0:
                aliq[i] = LiqAbsModel.liquid_water_absorption(denl[i], frq, tk[i])
            # compute ice absorption (db/km); convert non-zero value to np/km.
            if deni[i] > 0:
                aice[i] = np.dot(np.dot((8.18645 / wave), deni[i]), 0.000959553)
                aice[i] = np.dot(aice[i], db2np)

        return aliq, aice

    @staticmethod
    def clearsky_absorption(p=None, tk=None, e=None, frq=None, *args, **kwargs):
        """  Computes profiles of water vapor and dry air absorption for
        a given set of frequencies.  Subroutines H2O_xxx and O2_xxx
        contain the absorption model of Leibe and Layton [1987:
        Millimeter-wave properties of the atmosphere: laboratory studies
        and propagation modeling. NTIA Report 87-224, 74pp.] with oxygen
        interference coefficients from Rosenkranz [1988: Interference
        coefficients for overlapping oxygen lines in air.
        J. Quant. Spectrosc. Radiat. Transfer, 39, 287-97.]

        Args:
            p ([type], optional): pressure profile (mb). Defaults to None.
            tk ([type], optional): temperature profile (K). Defaults to None.
            e ([type], optional): vapor pressure profile (mb). Defaults to None.
            frq ([type], optional): frequency (GHz). Defaults to None.
            absmdl ([type], optional): Absorption model for WV (default 'ROS98')
            absmdl.wvres ([type], optional): wv resonant absorption
            absmdl.wvcnt ([type], optional): wv continuuum absorption

        Returns:
            (tuple):

                * awet [type]: water vapor absorption profile (np/km)
                * adry [type]: dry air absorption profile (np/km)

        See also:
            :py:meth:`absmodel.H2OAbsModel().h2o_rosen03_xxx`, :py:meth:`o2n2_rosen03_xxx`

        .. warning::
                * h2o_rosen03_xxx and o2n2_rosen03_xxx functions are missing!!!!!
                * rho, absmdl, absmdl.wvres and absmdl.wvcnt arguments are missing!!!!!!
        """

        nl = len(p)
        awet = np.zeros(p.shape)
        adry = np.zeros(p.shape)
        aO2 = np.zeros(p.shape)
        aN2 = np.zeros(p.shape)
        factor = np.dot(0.182, frq)
        db2np = np.dot(np.log(10.0), 0.1)
        for i in range(0, nl):
            # Compute inverse temperature parameter; convert wet and dry p to kpa.
            v = 300.0 / tk[i]
            ekpa = e[i] / 10.0
            pdrykpa = p[i] / 10.0 - ekpa

            # if H2OAbsModel.model == 'rose19sd':
            npp, ncpp = H2OAbsModel().h2o_rosen19_sd(pdrykpa, v, ekpa, frq)
            awet[i] = factor * (npp + ncpp) * db2np
            if O2AbsModel.model in ['rose19sd', 'rose19']:
                npp, ncpp = O2AbsModel().o2abs_rosen18(pdrykpa, v, ekpa, frq)
                aO2[i] = factor * npp * db2np
            if O2AbsModel.model in ['rose03', 'rose16']:
                npp, ncpp = O2AbsModel().o2abs_rosen19(pdrykpa, v, ekpa, frq)
                aO2[i] = factor * (npp + ncpp) * db2np
            # add N2 term
            if N2AbsModel.model not in ['rose03', 'rose16']:
                aN2[i] = N2AbsModel.n2_absorption(tk[i], np.dot(pdrykpa, 10), frq)

            if not N2AbsModel.model:
                raise ValueError('No model avalaible with this name: {} . Sorry...'.format('model'))

            adry[i] = aO2[i] + aN2[i]

        return awet, adry

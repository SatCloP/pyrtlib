# -*- coding: utf-8 -*-
"""
This class contains the absorption model used in pyrtlib.
"""

import numpy as np

from .linelist import cf, xcf, cs, xcs, fl, reftline, w0, w0s, xs, w2s, w2, sh, aair, xh, shs, aself, xhs, s1, b2, \
    reftcon, x
from .utils import dilec12, arange, dcerror


class AbsModel:
    """This class contains the absorption model used in pyrtlib.
    """
    model = ''
    """Model used to compute absorption"""

    @staticmethod
    def ab_liq(water=None, freq=None, temp=None, *args, **kwargs):
        """Computes Absorption In Nepers/Km By Suspended Water Droplets.

        Args:
            water ([type], optional): water in g/m3. Defaults to None.
            freq ([type], optional): frequency in GHz (Valid From 0 To 1000 Ghz). Defaults to None.
            temp ([type], optional): temperature in K. Defaults to None.

        Returns:
            [type]: [description]

        References
        ----------
        .. [1] Liebe, Hufford And Manabe, Int. J. Ir & Mm Waves V.12, Pp.659-675 (1991).
        .. [2] Liebe Et Al, Agard Conf. Proc. 542, May 1993.

        .. note::
            Revision history:

            * Pwr 8/3/92   Original Version
            * Pwr 12/14/98 Temp. Dependence Of Eps2 Eliminated To Agree With Mpm93 

        .. warning:: conversion of complex() function is missing
        """
        if water <= 0:
            abliq = 0
            return abliq

        if AbsModel.model == 'ros03':
            theta1 = 1.0 - 300.0 / temp
            eps0 = 77.66 - np.dot(103.3, theta1)
            eps1 = np.dot(0.0671, eps0)
            eps2 = 3.52
            fp = np.dot((np.dot(316.0, theta1) + 146.4), theta1) + 20.2
            fs = np.dot(39.8, fp)
            eps = (eps0 - eps1) / complex(1.0, freq / fp) + (eps1 - eps2) / complex(1.0, freq / fs) + eps2
        elif AbsModel.model == 'ros16':
            eps = dilec12(freq, temp)
        else:
            raise ValueError('No model avalaible with this name: {} . Sorry...'.format(AbsModel.model))

        re = (eps - 1.0) / (eps + 2.0)
        abliq = np.dot(np.dot(np.dot(- 0.06286, np.imag(re)), freq), water)

        return abliq

    @staticmethod
    def abs_N2(t=None, p=None, f=None, *args, **kwargs):
        """Collision-Induced Power Absorption Coefficient (Neper/Km) in air
        with modification of 1.34 to account for O2-O2 and O2-N2 collisions, as calculated by [3]

        5/22/02, 4/14/05 P.Rosenkranz

        Copyright (c) 2002 Massachusetts Institute of Technology

        Args:
            t ([type], optional): Temperature (K). Defaults to None.
            p ([type], optional): Dry Air Pressure (Mb). Defaults to None.
            f ([type], optional): Frequency (Ghz)(Valid 0-2000 Ghz). Defaults to None.

        Returns:
            [type]: [description]

        References
        ----------
        .. [1] See eq. 2.6 in Thermal Microwave Radiation - Applications for Remote Sensing (C. Maetzler, ed.) London, IET, 2006.
        .. [2] Borysow, A, and L. Frommhold, Astrophysical Journal, v.311, pp.1043-1057 (1986)
        .. [3] J.Boissoles, C.Boulet, R.H.Tipping, A.Brown, and Q.Ma, J. Quant. Spectros. Radiat. Trans. v.82, 505-516 (2003).
        """

        th = 300.0 / t
        fdepen = 0.5 + 0.5 / (1.0 + (f / 450.0) ** 2)
        if AbsModel.model in ['ros16', 'ros17', 'rose19sd']:
            l, m, n = 6.5e-14, 3.6, 1.34
        elif AbsModel.model == 'ros18':
            l, m, n = 9.9e-14, 3.22, 1
        elif AbsModel.model == 'ros03':
            l, m, n = 6.5e-14, 3.6, 1.29
        else:
            raise ValueError('No model avalaible with this name: {} . Sorry...'.format(AbsModel.model))

        bf = np.dot(np.dot(np.dot(np.dot(np.dot(np.dot(l, fdepen), p), p), f), f), th ** m)

        absN2 = np.dot(n, bf)

        return absN2

    @staticmethod
    def h2o_rosen19_sd(pdrykpa=None, vx=None, ekpa=None, frq=None, *args, **kwargs):
        """Compute absorption coef in atmosphere due to water vapor
        this version should not be used with a line list older than june 2018,
        nor the new list with an older version of this subroutine.
        Line parameters will be read from file h2o_list.asc; intensities
        should include the isotope abundance factors.
        this version uses a line-shape cutoff.

        Args:
            pdrykpa ([type], optional): [description]. Defaults to none.
            vx ([type], optional): [description]. Defaults to none.
            ekpa ([type], optional): [description]. Defaults to none.
            frq ([type], optional): [description]. Defaults to none.

        Returns:
            [type]: [description]

        References
        ----------
        .. [1] Rosenkranz, P.W.: Line-By-Line Microwave Radiative Transfer (Non-Scattering), Remote Sens. Code Library, Doi:10.21982/M81013, 2017
        """
        # nico: the best-fit voigt are given in koshelev et al. 2018, table 2 (rad,
        # mhz/torr). these correspond to w3(1) and ws(1) in h2o_list_r18 (mhz/mb)

        # read the list of parameters
        # h2o_sdlist_r19
        # h2o_sdlist_r20
        # this is the same as h2o_sdlist_r19 but the two coefficients w2air w2self at 22.2 ghz
        # (which were missing in h2o_sdlist_r19)

        # cyh ***********************************************
        db2np = np.dot(np.log(10.0), 0.1)
        rvap = np.dot(0.01, 8.31451) / 18.01528
        factor = np.dot(0.182, frq)
        t = 300.0 / vx
        p = np.dot((pdrykpa + ekpa), 10.0)
        rho = np.dot(ekpa, 10.0) / (np.dot(rvap, t))
        f = np.copy(frq)
        # cyh ***********************************************

        if rho <= 0.0:
            abh2o = 0.0
            npp = 0
            ncpp = 0

            return

        pvap = np.dot(rho, t) / 216.68

        pda = p - pvap
        den = np.dot(3.344e+16, rho)
        # continuum terms
        ti = reftcon / t
        # xcf and xcs include 3 for conv. to density & stimulated emission
        con = np.dot(
            np.dot(np.dot((np.dot(np.dot(cf, pda), ti ** xcf) + np.dot(np.dot(cs, pvap), ti ** xcs)), pvap), f), f)

        # nico 2019/03/18 *********************************************************
        # add resonances
        nlines = len(fl)
        ti = reftline / t
        tiln = np.log(ti)
        ti2 = np.exp(np.dot(2.5, tiln))

        sum = 0.0
        df = np.zeros((2, 1))
        for i in arange(1-1, nlines-1).reshape(-1):
            width0 = np.dot(np.dot(w0[i], pda), ti ** x[i]) + np.dot(np.dot(w0s[i], pvap), ti ** xs[i])
            width2 = np.dot(w2[i], pda) + np.dot(w2s[i], pvap)
            shiftf = np.dot(np.dot(np.dot(sh[i], pda), (1.0 - np.dot(aair[i], tiln))), ti ** xh[i])
            shifts = np.dot(np.dot(np.dot(shs[i], pvap), (1.0 - np.dot(aself[i], tiln))), ti ** xhs[i])
            shift = shiftf + shifts
            # nico: thus using the best-fit voigt (shift instead of shift0 and shift2)
            wsq = width0 ** 2
            s = np.dot(np.dot(s1[i], ti2), np.exp(np.dot(b2[i], (1.0 - ti))))
            df[0] = f - fl[i] - shift
            df[1] = f + fl[i] + shift
            base = width0 / (562500.0 + wsq)
            res = 0.0
            for j in arange(0, 1).reshape(-1):
                # if(i.eq.1 .and. j.eq.1 .and. abs(df(j)).lt.10.*width0) then
                # WIDTH2>0.0 & J==1 & abs(DF(J)) < 10*WIDTH0
                # width2 > 0 and j == 1 and abs(df[j]) < np.dot(10, width0)
                if width2 > 0 and j == 0 and abs(df[j]) < np.dot(10, width0):
                    # speed-dependent resonant shape factor
                    # double complex dcerror,xc,xrt,pxw,a
                    xc = complex((width0 - np.dot(1.5, width2)), df[j]) / width2
                    xrt = np.sqrt(xc)
                    pxw = np.dot(np.dot(1.77245385090551603, xrt), dcerror(-np.imag(xrt), np.real(xrt)))
                    sd = np.dot(2.0, (1.0 - pxw)) / width2
                    res = res + np.real(sd) - base
                else:
                    if abs(df[j]) < 750.0:
                        res = res + width0 / (df[j] ** 2 + wsq) - base
            sum = sum + np.dot(np.dot(s, res), (f / fl[i]) ** 2)
        # nico 2019/03/18 *********************************************************
        # cyh **************************************************************
        # separate the following original equ. into line and continuum
        # terms, and change the units from np/km to ppm
        # abh2o = .3183e-4*den*sum + con
        npp = (np.dot(np.dot(3.183e-05, den), sum) / db2np) / factor
        ncpp = (con / db2np) / factor
        # cyh *************************************************************

        return npp, ncpp

    @staticmethod
    def o2abs_rosen18_xxx(pdrykpa=None, vx=None, ekpa=None, frq=None, *args, **kwargs):
        """Returns power absorption coefficient due to oxygen in air,
        in nepers/km.  Multiply o2abs by 4.343 to convert to db/km.

        History:

        * 5/1/95  P. Rosenkranz
        * 11/5/97  P. Rosenkranz - 1- line modification.
        * 12/16/98 pwr - updated submm freq's and intensities from HITRAN96
        * 8/21/02  pwr - revised width at 425
        * 3/20/03  pwr - 1- line mixing and width revised
        * 9/29/04  pwr - new widths and mixing, using HITRAN intensities for all lines
        * 6/12/06  pwr - chg. T dependence of 1- line to 0.8
        * 10/14/08 pwr - moved isotope abundance back into intensities; added selected O16O18 lines.
        * 5/30/09  pwr - remove common block, add weak lines.
        * 12/18/14 pwr - adjust line broadening due to water vapor.

        Args:
            pdrykpa ([type], optional): [description]. Defaults to None.
            vx ([type], optional): [description]. Defaults to None.
            ekpa ([type], optional): [description]. Defaults to None.
            frq ([type], optional): [description]. Defaults to None.

        Returns:
            [type]: [description]

        References
        ----------
        .. [1] P.W. Rosenkranz, CHAP. 2 in ATMOSPHERIC REMOTE SENSING BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993) (http://hdl.handle.net/1721.1/68611).
        .. [3] G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, pp.282-287 (2003).
        .. [4] M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
        .. [5] M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
        .. [6] B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
        .. [7] D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
        .. [8] M.A. Koshelev et al, JQSRT, in press (2015). line intensities from HITRAN2004. non-resonant intensity from JPL catalog.

        .. note::
            * The mm line-width and mixing coefficients are from Tretyakov et al submm line-widths from Golubiatnikov & Krupnov (except 234 GHz from Drouin)
            * The same temperature dependence (X) is used for submillimeter line widths as in the 60 GHz band: (1/T)**X

        .. note::
            Nico: the only differences wrt to o2n2_rosen17_xxx are:

            * PRESWV = VAPDEN*TEMP/217. -> PRESWV = VAPDEN*TEMP/216.68 - this is now consistent with h2o
            * ABSN2_ros16(TEMP,PRES,FREQ) -> ABSN2_ros18(TEMP,PRESDA,FREQ) - since absn2.f takes in input P = DRY AIR PRESSURE (MB)
            * ABSN2 is now external
            * The continuum term is summed BEFORE O2ABS = max(O2ABS,0.)
        """

        # Nico 2016/11/30 *********************************************************
        # Here I imported the code I got from P. Rosenkranz on 2016/08/10, adapting
        # to RTE.

        # NB: NL=49
        # LINES ARE ARRANGED 1-,1+,...33-,33+ IN SPIN-ROTATION SPECTRUM;
        # BY FREQUENCY IN SUBMM SPECTRUM.
        # Nico F[GHz]
        f = np.array(
            [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.591, 59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
             56.9682, 62.4112, 56.3634, 62.998, 55.7838, 63.5685, 55.2214, 64.1278, 54.6712, 64.6789, 54.13, 65.2241,
             53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, 52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.431,
             50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.763, 487.2493, 566.8956, 715.3929, 731.1866, 773.8395,
             834.1455, 895.071])
        # Nico S(T_0)[Hz*cm2]
        s300 = np.array(
            [2.906e-15, 7.957e-16, 2.444e-15, 2.194e-15, 3.301e-15, 3.243e-15, 3.664e-15, 3.834e-15, 3.588e-15,
             3.947e-15, 3.179e-15, 3.661e-15, 2.59e-15, 3.111e-15, 1.954e-15, 2.443e-15, 1.373e-15, 1.784e-15,
             9.013e-16, 1.217e-15, 5.545e-16, 7.766e-16, 3.201e-16, 4.651e-16, 1.738e-16, 2.619e-16, 8.88e-17,
             1.387e-16, 4.272e-17, 6.923e-17, 1.939e-17, 3.255e-17, 8.301e-18, 1.445e-17, 3.356e-18, 6.049e-18,
             1.28e-18, 2.394e-18, 3.287e-17, 6.463e-16, 1.334e-17, 7.049e-15, 3.011e-15, 1.797e-17, 1.826e-15,
             2.193e-17, 1.153e-14, 3.974e-15, 2.512e-17])
        # Nico (Elow+hf)/kT_0 [unitless]
        be = np.array(
            [0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621, 0.91, 0.91, 1.255, 1.255, 1.654,
             1.654, 2.109, 2.109, 2.618, 2.618, 3.182, 3.182, 3.8, 3.8, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819,
             6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651, 0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136,
             0.141, 0.145, 0.201])

        # Nico gamma(T_0) [MHZ/mb == GHz/bar]
        wb300 = 0.56

        x = 0.8
        w300 = np.array(
            [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189,
             1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858,
             0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 1.65, 1.64, 1.64, 1.64, 1.6, 1.6, 1.6, 1.6, 1.62,
             1.47, 1.47])
        # nico y(t_0) [unitless]
        y300 = np.append(
            [- 0.036, 0.2547, - 0.3655, 0.5495, - 0.5696, 0.6181, - 0.4252, 0.3517, - 0.1496, 0.043, 0.064, - 0.1605,
             0.2906, - 0.373, 0.4169, - 0.4819, 0.4963, - 0.5481, 0.5512, - 0.5931, 0.6212, - 0.6558, 0.692, - 0.7208,
             0.7312, - 0.755, 0.7555, - 0.7751, 0.7914, - 0.8073, 0.8307, - 0.8431, 0.8676, - 0.8761, 0.9046, - 0.9092,
             0.9416, - 0.9423], np.tile(0.0, (1, 11)))
        # nico v(t_0) [unitless]
        v = np.append(
            [0.0079, - 0.0978, 0.0844, - 0.1273, 0.0699, - 0.0776, 0.2309, - 0.2825, 0.0436, - 0.0584, 0.6056, - 0.6619,
             0.6451, - 0.6759, 0.6547, - 0.6675, 0.6135, - 0.6139, 0.2952, - 0.2895, 0.2654, - 0.259, 0.375, - 0.368,
             0.5085, - 0.5002, 0.6206, - 0.6091, 0.6526, - 0.6393, 0.664, - 0.6475, 0.6729, - 0.6545, 0.68, - 0.66,
             0.685, - 0.665], np.tile(0.0, (1, 11)))
        # nico 2016/11/30 *********************************************************

        # cyh*** add the following lines *************************
        db2np = np.dot(np.log(10.0), 0.1)
        rvap = np.dot(0.01, 8.31451) / 18.01528
        factor = np.dot(0.182, frq)
        temp = 300.0 / vx
        pres = np.dot((pdrykpa + ekpa), 10.0)

        vapden = np.dot(ekpa, 10.0) / (np.dot(rvap, temp))
        freq = np.copy(frq)
        # cyh*****************************************************

        th = 300.0 / temp
        th1 = th - 1.0
        b = th ** x
        preswv = np.dot(vapden, temp) / 216.68

        presda = pres - preswv

        den = np.dot(0.001, (np.dot(presda, b) + np.dot(np.dot(1.2, preswv), th)))

        dfnr = np.dot(wb300, den)

        # nico intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
        sum = np.dot(np.dot(np.dot(1.584e-17, freq), freq), dfnr) / (
            np.dot(th, (np.dot(freq, freq) + np.dot(dfnr, dfnr))))
        # cyh **************************************************************

        nlines = len(f)
        for k in arange(1-1, nlines-1).reshape(-1):
            df = np.dot(w300[k], den)
            fcen = f[k]
            y = np.dot(den, (y300[k] + np.dot(v[k], th1)))
            strr = np.dot(s300[k], np.exp(np.dot(-be[k], th1)))
            sf1 = (df + np.dot((freq - fcen), y)) / ((freq - fcen) ** 2 + np.dot(df, df))
            sf2 = (df - np.dot((freq + fcen), y)) / ((freq + fcen) ** 2 + np.dot(df, df))
            sum = sum + np.dot(np.dot(strr, (sf1 + sf2)), (freq / f[k]) ** 2)

        # o2abs = .5034e12*sum*presda*th^3/3.14159;
        # .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
        # nico the following computes n/pi*sum*th^2 (see notes)
        # nico n/pi = a/(pi*k*t_0) * pda * th
        # nico a/(pi*k*t_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont
        # for units conversion (pa->hpa, hz->ghz, m->km)
        # nico pa2hpa=1e-2; hz2ghz=1e-9; m2cm=1e2; m2km=1e-3; pa2hpa^-1 * hz2ghz * m2cm^-2 * m2km^-1 = 1e-8
        # nico th^3 = th(from ideal gas law 2.13) * th(from the mw approx of stimulated emission 2.16 vs. 2.14) *
        # th(from the partition sum 2.20)
        o2abs = np.dot(np.dot(np.dot(1.6097e+11, sum), presda), th ** 3)
        o2abs = np.maximum(o2abs, 0.0)
        # cyh *** ********************************************************
        # separate the equ. into line and continuum
        # terms, and change the units from np/km to ppm

        npp = np.copy(o2abs)

        # nico intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
        ncpp = np.dot(np.dot(np.dot(1.584e-17, freq), freq), dfnr) / (
            np.dot(th, (np.dot(freq, freq) + np.dot(dfnr, dfnr))))
        # 20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
        # nico a/(pi*k*t_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont for units conversion (pa->hpa, hz->ghz)
        # nico pa2hpa=1e-2; hz2ghz=1e-9; m2cm=1e2; m2km=1e-3; pa2hpa^-1 * hz2ghz * m2cm^-2 * m2km^-1 = 1e-8
        # nico th^3 = th(from ideal gas law 2.13) * th(from the mw approx of stimulated emission 2.16 vs. 2.14) * th(from the partition sum 2.20)
        ncpp = np.dot(np.dot(np.dot(1.6097e+11, ncpp), presda), th ** 3)

        # change the units from np/km to ppm
        npp = (npp / db2np) / factor

        ncpp = (ncpp / db2np) / factor
        # cyh ************************************************************

        return npp, ncpp

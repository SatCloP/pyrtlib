# -*- coding: utf-8 -*-
"""
This class contains the absorption model used in pyrtlib.
"""

import numpy as np

from .linelist import h2o_linelist as h20ll
from .linelist import o2_linelist
from .utils import dilec12, dcerror


class AbsModel(object):
    """This is an abstraction class to define the absorption model.
    """

    def __init__(self, model):
        self._model = model
        """Model used to compute absorption"""

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        if model and isinstance(model, str):
            self._model = model
        else:
            raise ValueError("Please enter a valid absorption model")


class LiqAbsModel(AbsModel):
    """This class contains the absorption model used in pyrtlib.
    """

    @staticmethod
    def liquid_water_absorption(water=None, freq=None, temp=None, *args, **kwargs):
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

        if LiqAbsModel.model == 'ros03':
            theta1 = 1.0 - 300.0 / temp
            eps0 = 77.66 - np.dot(103.3, theta1)
            eps1 = np.dot(0.0671, eps0)
            eps2 = 3.52
            fp = np.dot((np.dot(316.0, theta1) + 146.4), theta1) + 20.2
            fs = np.dot(39.8, fp)
            eps = (eps0 - eps1) / complex(1.0, freq / fp) + (eps1 - eps2) / complex(1.0, freq / fs) + eps2
        elif LiqAbsModel.model == 'ros16':
            eps = dilec12(freq, temp)
        else:
            raise ValueError('[AbsLiq] No model available with this name: {} . Sorry...'.format(LiqAbsModel.model))

        re = (eps - 1.0) / (eps + 2.0)

        abliq = -0.06286 * np.imag(re) * freq * water

        return abliq


class N2AbsModel(AbsModel):
    """This class contains the absorption model used in pyrtlib.
    """

    @staticmethod
    def n2_absorption(t=None, p=None, f=None, *args, **kwargs):
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
        if N2AbsModel.model in ['ros16', 'ros17', 'rose19sd']:
            l, m, n = 6.5e-14, 3.6, 1.34
        elif N2AbsModel.model == 'ros18':
            l, m, n = 9.9e-14, 3.22, 1
        elif N2AbsModel.model == 'ros03':
            l, m, n = 6.5e-14, 3.6, 1.29
        else:
            raise ValueError('[AbsN2] No model available with this name: {} . Sorry...'.format(N2AbsModel.model))

        bf = l * fdepen * p * p * f * f * th ** m

        absN2 = n * bf

        return absN2


class H2OAbsModel(AbsModel):
    """This class contains the absorption model used in pyrtlib.
    """

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
        db2np = np.log(10.0) * 0.1
        rvap = (0.01 * 8.31451) / 18.01528
        factor = 0.182 * frq
        t = 300.0 / vx
        p = (pdrykpa + ekpa) * 10.0
        rho = ekpa * 10.0 / (rvap * t)
        f = np.copy(frq)
        # cyh ***********************************************

        if rho <= 0.0:
            abh2o = 0.0
            npp = 0
            ncpp = 0

            return

        pvap = (rho * t) / 216.68

        pda = p - pvap
        den = 3.344e+16 * rho
        # continuum terms
        ti = h20ll.reftcon / t
        # xcf and xcs include 3 for conv. to density & stimulated emission

        con = (h20ll.cf * pda * ti ** h20ll.xcf + h20ll.cs * pvap * ti ** h20ll.xcs) * pvap * f * f

        # nico 2019/03/18 *********************************************************
        # add resonances
        nlines = len(h20ll.fl)
        ti = h20ll.reftline / t
        tiln = np.log(ti)
        ti2 = np.exp(2.5 * tiln)

        sum = 0.0
        df = np.zeros((2, 1))
        for i in range(0, nlines):
            width0 = h20ll.w0[i] * pda * ti ** h20ll.x[i] + h20ll.w0s[i] * pvap * ti ** h20ll.xs[i]
            width2 = h20ll.w2[i] * pda + h20ll.w2s[i] * pvap

            shiftf = h20ll.sh[i] * pda * (1. - h20ll.aair[i] * tiln) * ti ** h20ll.xh[i]
            shifts = h20ll.shs[i] * pvap * (1. - h20ll.aself[i] * tiln) * ti ** h20ll.xhs[i]
            shift = shiftf + shifts
            # nico: thus using the best-fit voigt (shift instead of shift0 and shift2)
            wsq = width0 ** 2
            s = h20ll.s1[i] * ti2 * np.exp(h20ll.b2[i] * (1. - ti))
            df[0] = f - h20ll.fl[i] - shift
            df[1] = f + h20ll.fl[i] + shift
            base = width0 / (562500.0 + wsq)
            res = 0.0
            for j in range(0, 1):
                # if(i.eq.1 .and. j.eq.1 .and. abs(df(j)).lt.10.*width0) then
                # WIDTH2>0.0 & J==1 & abs(DF(J)) < 10*WIDTH0
                # width2 > 0 and j == 1 and abs(df[j]) < np.dot(10, width0)
                if width2 > 0 and j == 0 and abs(df[j]) < np.dot(10, width0):
                    # speed-dependent resonant shape factor
                    # double complex dcerror,xc,xrt,pxw,a
                    xc = complex((width0 - np.dot(1.5, width2)), df[j]) / width2
                    xrt = np.sqrt(xc)

                    pxw = 1.77245385090551603 * xrt * dcerror(-np.imag(xrt), np.real(xrt))
                    sd = 2.0 * (1.0 - pxw) / width2
                    res = res + np.real(sd) - base
                else:
                    if abs(df[j]) < 750.0:
                        res = res + width0 / (df[j] ** 2 + wsq) - base

            sum = sum + s * res * (f / h20ll.fl[i]) ** 2
        # nico 2019/03/18 *********************************************************
        # cyh **************************************************************
        # separate the following original equ. into line and continuum
        # terms, and change the units from np/km to ppm
        # abh2o = .3183e-4*den*sum + con

        npp = (3.183e-05 * den * sum / db2np) / factor
        ncpp = (con / db2np) / factor
        # cyh *************************************************************

        return npp, ncpp


class O2AbsModel(AbsModel):
    """This class contains the absorption model used in pyrtlib.
    """

    @staticmethod
    def o2abs_rosen18(pdrykpa=None, vx=None, ekpa=None, frq=None, *args, **kwargs):
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

        o2ll = o2_linelist(O2AbsModel.model)

        # cyh*** add the following lines *************************
        db2np = np.log(10.0) * 0.1
        rvap = (0.01 * 8.31451) / 18.01528
        factor = 0.182 * frq
        temp = 300.0 / vx
        pres = (pdrykpa + ekpa) * 10.0

        vapden = ekpa * 10.0 / (rvap * temp)
        freq = np.copy(frq)
        # cyh*****************************************************

        th = 300.0 / temp
        th1 = th - 1.0
        b = th ** o2ll.x
        preswv = (vapden * temp) / 216.68
        presda = pres - preswv
        den = 0.001 * (presda * b + 1.2 * preswv * th)
        dfnr = o2ll.wb300 * den

        # nico intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
        # 1.584e-17*freq*freq*dfnr/(th*(freq*freq + dfnr*dfnr))
        sum = 1.584e-17 * freq * freq * dfnr / (th * (freq * freq + dfnr * dfnr))
        # cyh **************************************************************

        nlines = len(o2ll.f)
        for k in range(0, nlines):
            df = o2ll.w300[k] * den
            fcen = o2ll.f[k]

            y = den * (o2ll.y300[k] + o2ll.v[k] * th1)
            strr = o2ll.s300[k] * np.exp(-o2ll.be[k] * th1)
            sf1 = (df + (freq - fcen) * y) / ((freq - fcen) ** 2 + df * df)
            sf2 = (df - (freq + fcen) * y) / ((freq + fcen) ** 2 + df * df)
            sum = sum + strr * (sf1 + sf2) * (freq / o2ll.f[k]) ** 2

        # o2abs = .5034e12*sum*presda*th^3/3.14159;
        # .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
        # nico the following computes n/pi*sum*th^2 (see notes)
        # nico n/pi = a/(pi*k*t_0) * pda * th
        # nico a/(pi*k*t_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont
        # for units conversion (pa->hpa, hz->ghz, m->km)
        # nico pa2hpa=1e-2; hz2ghz=1e-9; m2cm=1e2; m2km=1e-3; pa2hpa^-1 * hz2ghz * m2cm^-2 * m2km^-1 = 1e-8
        # nico th^3 = th(from ideal gas law 2.13) * th(from the mw approx of stimulated emission 2.16 vs. 2.14) *
        # th(from the partition sum 2.20)

        o2abs = 1.6097e+11 * sum * presda * th ** 3
        o2abs = np.maximum(o2abs, 0.0)
        # cyh *** ********************************************************
        # separate the equ. into line and continuum
        # terms, and change the units from np/km to ppm

        npp = np.copy(o2abs)

        # nico intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17

        ncpp = 1.584e-17 * freq * freq * dfnr / (th * (freq * freq + dfnr * dfnr))
        # 20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
        # nico a/(pi*k*t_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont
        # for units conversion (pa->hpa, hz->ghz)
        # nico pa2hpa=1e-2; hz2ghz=1e-9; m2cm=1e2; m2km=1e-3; pa2hpa^-1 * hz2ghz * m2cm^-2 * m2km^-1 = 1e-8
        # nico th^3 = th(from ideal gas law 2.13) * th(from the mw approx of stimulated emission 2.16 vs. 2.14)
        # * th(from the partition sum 2.20)
        ncpp = 1.6097e+11 * ncpp * presda * th ** 3

        # change the units from np/km to ppm
        npp = (npp / db2np) / factor

        ncpp = (ncpp / db2np) / factor
        # cyh ************************************************************

        return npp, ncpp

    @staticmethod
    def o2abs_rosen19(pdrykpa=None, vx=None, ekpa=None, frq=None, *args, **kwargs):
        """Returns power absorption coefficient due to oxygen in air,
        in nepers/km. multiply o2abs2 by 4.343 to convert to db/km.

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
            pdrykpa ([type], optional): [description]. Defaults to None.
            vx ([type], optional): [description]. Defaults to None.
            ekpa ([type], optional): [description]. Defaults to None.
            frq ([type], optional): [description]. Defaults to None.

        Returns:
            [type]: [description]

        References
        ----------
        .. [1] P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993)(http://hdl.handle.net/1721.1/68611).
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
            line widths as in the 60 GHz band: (1/T)**X (Koshelev et al 2016)
         3. The sign of DNU in the shape factor is corrected.
        """

        o2ll = o2_linelist(O2AbsModel.model)

        # *** add the following lines *************************
        db2np = np.log(10.0) * 0.1
        rvap = (0.01 * 8.31451) / 18.01528
        factor = 0.182 * frq
        temp = 300.0 / vx
        pres = (pdrykpa + ekpa) * 10.0
        vapden = ekpa * 10.0 / (rvap * temp)
        freq = np.copy(frq)
        # *****************************************************

        th = 300.0 / temp
        th1 = th - 1.0
        b = th ** o2ll.x
        preswv = vapden * temp / 216.68
        presda = pres - preswv
        den = 0.001 * (presda * b + 1.2 * preswv * th)
        dfnr = o2ll.wb300 * den
        pe2 = den * den
        # nico intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17

        sum = 1.584e-17 * freq * freq * dfnr / (th * (freq * freq + dfnr * dfnr))
        nlines = len(o2ll.f)
        for k in range(0, nlines):
            y = den * (o2ll.y0[k] + o2ll.y1[k] * th1)
            dnu = pe2 * (o2ll.dnu0[k] + o2ll.dnu1[k] * th1)
            gfac = 1. + pe2 * (o2ll.g0[k] + o2ll.g1[k] * th1)
            df = o2ll.w300[k] * den
            strr = o2ll.s300[k] * np.exp(-o2ll.be[k] * th1)
            del1 = freq - o2ll.f[k] - dnu
            del2 = freq + o2ll.f[k] + dnu
            d1 = del1 * del1 + df * df
            d2 = del2 * del2 + df * df
            sf1 = (df * gfac + del1 * y) / d1
            sf2 = (df * gfac - del2 * y) / d2
            sum = sum + strr * (sf1 + sf2) * (freq / o2ll.f[k]) ** 2

        # .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11

        o2abs = 1.6097e+11 * sum * presda * th ** 3
        o2abs = 1.004 * np.maximum(o2abs, 0.0)

        # *** ********************************************************
        # separate the equ. into line and continuum
        # terms, and change the units from np/km to ppm

        npp = np.copy(o2abs)

        # nico intensities of the non-resonant transitions for o16-o16 and o16-o18, from jpl's line compilation
        # 1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17

        ncpp = 1.584e-17 * freq * freq * dfnr / (th * (freq * freq + dfnr * dfnr));
        #  .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
        # nico a/(pi*k*t_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont
        # for units conversion (pa->hpa, hz->ghz)
        # nico pa2hpa=1e-2; hz2ghz=1e-9; m2cm=1e2; m2km=1e-3; pa2hpa^-1 * hz2ghz * m2cm^-2 * m2km^-1 = 1e-8
        # nico th^3 = th(from ideal gas law 2.13) * th(from the mw approx of stimulated emission 2.16 vs. 2.14) *
        # th(from the partition sum 2.20)
        ncpp = 1.6097e11 * ncpp * presda * th ** 3  # nico: n/pi*sum0

        # change the units from np/km to ppm
        npp = (npp / db2np) / factor

        ncpp = (ncpp / db2np) / factor
        # ************************************************************

        return npp, ncpp

# -*- coding: utf-8 -*-
"""
This class contains the absorption model used in pyrtlib.
"""

import numpy as np
from .utils import dilec12


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
        .. [1] Liebe, Hufford And Manabe, Int. J. Ir & Mm Waves V.12, Pp.659-675 (1991);  
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
        """Absn2 = Collision-Induced Power Absorption Coefficient (Neper/Km) in air
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
        .. [1] See eq. 2.6 in Thermal Microwave Radiation - Applications
                for Remote Sensing (C. Maetzler, ed.) London, IET, 2006.
        .. [2] Borysow, A, and L. Frommhold,
                Astrophysical Journal, v.311, pp.1043-1057 (1986)
        .. [3] J.Boissoles, C.Boulet, R.H.Tipping, A.Brown, and Q.Ma,
                J. Quant. Spectros. Radiat. Trans. v.82, 505-516 (2003).
        """

        th = 300.0 / t
        fdepen = 0.5 + 0.5 / (1.0 + (f / 450.0) ** 2)
        if AbsModel.model == 'ros16':
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

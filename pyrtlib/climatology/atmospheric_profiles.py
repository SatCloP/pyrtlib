# -*- coding: utf-8 -*-
"""
This calss contains the AFGL Atmospheric Constituent Profiles (0-120km).
"""
import os
from typing import Tuple, Dict
import numpy as np


class AtmosphericProfiles:
    """AFGL Atmospheric Constituent Profiles (0-120km)

    Each of these profile contains data at 50 atmospheric levels.
    Altitude (km), Pressure (mb), Density (cm-3), Molec. densities (ppmv):

    * 0 (H2O),
    * 1 (CO2),
    * 2 (O3),
    * 3 (N2O),
    * 4 (CO),
    * 5 (CH4),
    * 6 (O2)

    Plus suplimental profiles where available.
    The last set of data sets are constituent profiles of molecular
    densities (ppmv) for the minor absorbing atmospheric gases.

    References
    ----------

    .. [1] [ANDERSON]_

    Examples:
        .. code-block:: python

            >>> from pyrtlib.climatology import AtmosphericProfiles as atmp
            >>> atmp.atm_profiles()
            {0: 'Tropical',
            1: 'Midlatitude Summer',
            2: 'Midlatitude Winter',
            3: 'Subarctic Summer',
            4: 'Subarctic Winter',
            5: 'US Standard'}
            >>> atmp.TROPICAL, atmp.H2O
            (0, 0)
    """

    TROPICAL = 0
    MIDLATITUDE_SUMMER = 1
    MIDLATITUDE_WINTER = 2
    SUBARCTIC_SUMMER = 3
    SUBARCTIC_WINTER = 4
    US_STANDARD = 5

    H2O = 0
    CO2 = 1
    O3 = 2
    N2O = 3
    CO = 4
    CH4 = 5
    O2 = 6
    # MINOR
    NO = 7
    SO2 = 8
    NO2 = 9
    NH3 = 10
    HNO3 = 11
    OH = 12
    HF = 13
    HCL = 14
    HBR = 15
    HI = 16
    CLO = 17
    OCS = 18
    H2CO = 19
    HOCL = 20
    N2 = 21
    HCN = 22
    CH3CL = 23
    H2O2 = 24
    C2H2 = 25
    C2H6 = 26
    PH3 = 27
    # TRACE
    COF2 = 28
    SF6 = 29
    H2S = 30
    CCL3F = 31
    CCL2F2 = 32
    CCLF3 = 33
    CF4 = 34
    CHCl2F = 35
    CHCLF2 = 36
    C2CL3F3 = 37
    C2CL2F4 = 38
    C2CLF5 = 39
    CCL4 = 40
    CLONO2 = 41
    N2O5 = 42
    HNO4 = 43
    AIR = 99

    @staticmethod
    def atm_profiles() -> Dict[int, str]:
        """Convenient function to ger the list of the buolt-in  Atmospheric profiles.

        Returns:
            Dict[int, str]: A dictionary of standard profiles atmospheric.
        """
        return {AtmosphericProfiles.TROPICAL: 'Tropical',
                AtmosphericProfiles.MIDLATITUDE_SUMMER: 'Midlatitude Summer',
                AtmosphericProfiles.MIDLATITUDE_WINTER: 'Midlatitude Winter',
                AtmosphericProfiles.SUBARCTIC_SUMMER: 'Subarctic Summer',
                AtmosphericProfiles.SUBARCTIC_WINTER: 'Subarctic Winter',
                AtmosphericProfiles.US_STANDARD: 'US Standard'}

    @staticmethod
    def gl_atm(atm: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Returns the Atmopshere profile.

        This method contains 6 model profiles:

        +--------+---------------------+
        | option | model               |
        +--------+---------------------+
        | 1      |  Tropical           |
        +--------+---------------------+
        | 2      |  Midlatitude Summer |
        +--------+---------------------+
        | 3      |  Midlatitude Winter |
        +--------+---------------------+
        | 4      |  Subarctic Summer   |
        +--------+---------------------+
        | 5      |  Subarctic Winter   |
        +--------+---------------------+
        | 6      |  U.S. Standard      |
        +--------+---------------------+

        Args:
            option (int): the atmosphere profile.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]: 
            * a (numpy.ndarray): Altitudes (km) (50x1)
            * p (numpy.ndarray): Pressure (mbar)
            * d (numpy.ndarray): Total density (cm-3)
            * t (numpy.ndarray): Temperature (K)
            * md (numpy.ndarray): Molecular densities (ppmv)

        Example:
            >>> from pyrtlib.climatology import AtmosphericProfiles as atmp
            >>> z, p, _, t, md = atmp.gl_atm(atmp.US_STANDARD)
            >>> md[:, atmp.H2O]
            array([7.745e+03, 6.071e+03, 4.631e+03, 3.182e+03, 2.158e+03, 1.397e+03,
                  9.254e+02, 5.720e+02, 3.667e+02, 1.583e+02, 6.996e+01, 3.613e+01,
                  1.906e+01, 1.085e+01, 5.927e+00, 5.000e+00, 3.950e+00, 3.850e+00,
                  3.825e+00, 3.850e+00, 3.900e+00, 3.975e+00, 4.065e+00, 4.200e+00,
                  4.300e+00, 4.425e+00, 4.575e+00, 4.725e+00, 4.825e+00, 4.900e+00,
                  4.950e+00, 5.025e+00, 5.150e+00, 5.225e+00, 5.250e+00, 5.225e+00,
                  5.100e+00, 4.750e+00, 4.200e+00, 3.500e+00, 2.825e+00, 2.050e+00,
                  1.330e+00, 8.500e-01, 5.400e-01, 4.000e-01, 3.400e-01, 2.800e-01,
                  2.400e-01, 2.000e-01])

        .. note:: adapted from glatm.dat.  DCT 3/26/97
        """

        path = os.path.dirname(os.path.abspath(__file__))
        option = atm + 1
        # MODEL 1.  TROPICAL
        if option == 1:
            #
            # Latitude (deg)
            # 15.0
            mtx = np.loadtxt(os.path.join(path, "tropical.dat"))

        # MODEL 2. MIDLATITUDE SUMMER
        elif option == 2:
            #
            # Latitude (deg)
            # 45.0
            mtx = np.loadtxt(os.path.join(path, "midlatitude_summer.dat"))

        # MODEL 3. MIDLATITUDE WINTER
        elif option == 3:
            #
            # Latitude (deg)
            # 45.0
            mtx = np.loadtxt(os.path.join(path, "midlatitude_winter.dat"))

        # MODEL 4. SUBARCTIC SUMMER
        elif option == 4:
            #
            # Latitude (deg)
            # 60.0
            #
            mtx = np.loadtxt(os.path.join(path, "subarctic_summer.dat"))

        # MODEL 5. SUBARCTIC WINTER
        elif option == 5:
            #
            # Latitude (deg)
            # 60.0
            mtx = np.loadtxt(os.path.join(path, "subarctic_winter.dat"))

        # MODEL 6. U.S. STANDARD
        elif option == 6:
            #
            # Latitude (deg)
            # 45.5397
            mtx = np.loadtxt(os.path.join(path, "us_standard.dat"))

        # Altitude (km)
        a = mtx[:, 0]
        # Pressue (mb)
        p = mtx[:, 1]
        # Temperature (K)
        t = mtx[:, 3]
        # Density (cm-3)
        d = mtx[:, 2]
        # 1 :H2O (ppmv)
        g1 = mtx[:, 4]
        # 2 :CO2 (ppmv)
        g2 = mtx[:, 5]
        # 3 :O3 (ppmv)
        g3 = mtx[:, 6]
        # 4 : N2O (ppmv)
        g4 = mtx[:, 7]
        # 5 :CO (ppmv)
        g5 = mtx[:, 8]
        # 6 :CH4 (ppmv)
        g6 = mtx[:, 9]
        # 7 :O2 (ppmv)
        g7 = mtx[:, 10]

        md = np.column_stack((g1, g2, g3, g4, g5, g6, g7))

        return a, p, d, t, md

    @staticmethod
    def gl_atm_minor(gas_minor: int) -> np.ndarray:
        """Returns the minor gas profiles (gas ID's 8-28)

        Args:
            gas_minor (int): HITRAN gas ID #'s (#gases x 1)

        Returns:
            numpy.ndarray: molecular densities (ppmv) (50x #gases)

        Example:

            >>> from pyrtlib.climatology import AtmosphericProfiles as atmp
            >>> atmp.gl_atm_minor(atmp.NO)
            array([3.00e-04, 3.00e-04, 3.00e-04, 3.00e-04, 3.00e-04, 3.00e-04,
                   3.00e-04, 3.00e-04, 3.00e-04, 3.00e-04, 3.00e-04, 3.00e-04,
                   3.00e-04, 2.99e-04, 2.95e-04, 2.83e-04, 2.68e-04, 2.52e-04,
                   2.40e-04, 2.44e-04, 2.55e-04, 2.77e-04, 3.07e-04, 3.60e-04,
                   4.51e-04, 6.85e-04, 1.28e-03, 2.45e-03, 4.53e-03, 7.14e-03,
                   9.34e-03, 1.12e-02, 1.19e-02, 1.17e-02, 1.10e-02, 1.03e-02,
                   1.01e-02, 1.01e-02, 1.03e-02, 1.15e-02, 1.61e-02, 2.68e-02,
                   7.01e-02, 2.13e-01, 7.12e-01, 2.08e+00, 4.50e+00, 7.98e+00,
                   1.00e+01, 1.00e+01])
        """

        path = os.path.dirname(os.path.abspath(__file__))
        mtx = np.loadtxt(os.path.join(path, "gas_minor.dat"))
        gas_minor -= 7
        # CONSTITUENT PROFILES FOR THE MINOR ABSORBING ATMOSPHERIC GASES
        #
        # MINGAS
        g = mtx[:, gas_minor]

        return g

    @staticmethod
    def gl_atm_trace(gas_trace: int) -> np.ndarray:
        """Returns the trace gas profiles (ID's 29-31,51-63)

        Args:
            gas_trace (int): HITRAN gas ID

        Returns:
            numpy.ndarray: molecular densities (ppmv) (50x #gases)

        Example:

            >>> from pyrtlib.climatology import AtmosphericProfiles as atmp
            >>> atmp.gl_atm_trace(atmp.H2S)
            array([1.00e-04, 5.48e-05, 3.00e-05, 2.45e-05, 2.00e-05, 1.61e-05,
                  1.30e-05, 1.14e-05, 1.00e-05, 8.37e-06, 7.00e-06, 4.58e-06,
                  3.00e-06, 1.73e-06, 1.00e-06, 5.48e-07, 3.00e-07, 1.73e-07,
                  1.00e-07, 5.48e-08, 3.00e-08, 1.73e-08, 1.00e-08, 5.48e-09,
                  3.00e-09, 1.73e-09, 1.00e-09, 1.00e-10, 1.00e-11, 1.00e-12,
                  1.00e-13, 1.00e-14, 1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15,
                  1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15,
                  1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15, 1.00e-15,
                  1.00e-15, 1.00e-15])

        References
        ----------

            .. [1] Smith, M.A.H., Compilation of atmospheric gas concentration profiles from 0 to 50km, NASA Tech.Mem. 83289, 1982
            .. [2] Rinsland, C.P., et al., JGR 94 D15 18341-18349, 1989
            .. [3] Rinsland, C.P., et al., JGR 95 D10 16477-16490, 1990

        """

        # ************************************************************************
        # ATMOSPHERIC PROFILES OF THE MINOR GASES (ID>28)
        # REFERENCES
        # ATMOS with no further reference refers to the profile
        # set distributed by diskette
        # Smith, M.A.H., Compilation of atmospheric gas
        # concentration profiles from 0 to 50km,
        # NASA Tech.Mem. 83289, 1982
        # Rinsland, C.P., et al., JGR 94 D15 18341-18349, 1989
        #
        # Rinsland, C.P., et al., JGR 95 D10 16477-16490, 1990
        # Copy of FASCODE XMLATM BLOCK DATA
        # Data from Mark Allen at JPL/private communication to Gail Anderson
        # ************************************************************************

        path = os.path.dirname(os.path.abspath(__file__))
        mtx = np.loadtxt(os.path.join(path, "gas_trace.dat"))
        gas_trace -= 28
        #  COF2 & SF6
        #   0 - 14   ATMOS SS Rinsland 1990: ConsAt 14 value
        #  14 -120   ATMOS SS Rinsland 1990: LinExLogMix until 1.0E-15

        #  H2S, CCLF3  (CFC13), CHCl2F (CFC21), C2CLF5  (CFC115)
        #   0 - 120   Smith 1982: until 1.0E-15
        #  CCL3F (F11), CCL2F2 (CFC12), CHCLF2 (CFC22), C2CL3F3 (CFC113), C2CL2F4 (CFC114), CCL4, CLONO2, N2O5, HNO4
        #   0 -120   Mark Allen/FASCODE
        #  CF4 (CFC14)
        #   0 - 10   ATMOS SS: ConsAt 10.5 value
        #  11 - 30   ATMOS SS: LinExLogMix
        #  30 - 75   ATMOS SS: ConsAt 31.5 value
        #  75 -120   1.0E-15
        g = mtx[:, gas_trace]

        return g

"""
This module contains the utils functions.
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

from typing import Tuple, Optional

import numpy as np


def import_lineshape(name):
    """ Import a named object from a module in the context of this function.
    """
    try:
        module = __import__('pyrtlib.lineshape', globals(), locals(), [name])
    except ImportError:
        return None
    return vars(module)[name]


def constants(string: str) -> Tuple[float, str]:
    """
    This routine will provide values and units for all the
    universal constants that I needed in my work.

    Nico, 2000

    History
    2003/06/06 - Added 'Rdry' and 'Rwatvap'
    2005/12/10 - Added 'Tcosmicbkg'
    2021/03/09 - Added 'avogadro', 'gravity'

    Args:
        string (str): String specifying which constant is needed. Defaults to None.

    +---------------+-------------------------------------------+
    | string        | description                               |
    +---------------+-------------------------------------------+
    | 'avogadro'    |  Avogadro number [mol-1]                  |
    +---------------+-------------------------------------------+
    | 'boltzmann'   |  Boltzmann constant [J K-1]               |
    +---------------+-------------------------------------------+
    | 'EarthRadius' |  Earth radius [km]                        |
    +---------------+-------------------------------------------+
    | 'light'       |  Light speed [m s-1]                      |
    +---------------+-------------------------------------------+
    | 'Np2dB'       |  Neper to Decibel [dB/Np]                 |
    +---------------+-------------------------------------------+
    | 'planck'      |  Planck constant [J Hz-1]                 |
    +---------------+-------------------------------------------+
    | 'Rdry'        |  Gas constant of dry air [J kg-1 K-1]     |
    +---------------+-------------------------------------------+
    | 'Rwatvap'     |  Gas constant of water vapor [J kg-1 K-1] |
    +---------------+-------------------------------------------+
    | 'R'           |  Gas constant [J mol-1 K-1]               |
    +---------------+-------------------------------------------+
    | 'Tcosmicbkg'  |  Cosmic Background Temperature [K]        |
    +---------------+-------------------------------------------+

    Raises:
        ValueError: [description]

    Returns:
        (tuple) : Numerical Value of the asked constant and string specifying which units are used

    References
    ----------
    .. [1] P.J. Mohr, B.N. Taylor, and D.B. Newell (2015), "The 2014 CODATA Recommended
            Values of the Fundamental Physical Constants" (Web Version 7.0), http://physics.nist.gov/cuu/index.html
            Values as of 11/12/2015
    """

    if string == 'avogadro':
        nA = 6.022140857e+23
        units = '[mol-1]'
        out = np.copy(nA)
    elif string == 'boltzmann':
        K = np.dot(1.380658, 1e-23)
        units = '[J K-1]'
        out = np.copy(K)
    elif string == 'EarthRadius':
        R = 6370.949
        units = '[Km]'
        out = np.copy(R)
    elif string == 'gravity':
        g = 9.80665
        units = '[m s-2]'
        out = np.copy(g)
    elif string == 'light':
        c = 299792458
        units = '[m s-1]'
        out = np.copy(c)
    elif string == 'Np2dB':
        Np2dB = np.dot(10, np.log10(np.exp(1)))
        units = '[dB/Np]'
        out = np.copy(Np2dB)
    elif string == 'planck':
        h = np.dot(6.6260755, 1e-34)
        units = '[J Hz-1]'
        out = np.copy(h)
    elif string == 'Rdry':
        Rd = 287.04
        units = '[J kg-1 K-1]'
        out = np.copy(Rd)
    elif string == 'Rwatvap':
        Rv = 461.5
        units = '[J kg-1 K-1]'
        out = np.copy(Rv)
    elif string == 'Tcosmicbkg':
        # Tcos = 2.736; # +/- 0.017 [K] Cosmic Background Temperature,
        # from Janssen, Atmospheric Remote Sensing by Microwave Radiometry, pag.12
        Tcos = 2.728
        units = '[K]'
        out = np.copy(Tcos)
    elif string == 'R':
        Rgas = 8.31446261815324
        units = '[J mol-1 K-1]'
        out = np.copy(Rgas)
    else:
        raise ValueError(
            'No constant avalaible with this name: {} . Sorry...'.format(
                string))

    return out, units


def gas_mass(gasid: int) -> float:
    """
    Returns the mass of the HITRAN gas ID
    DCT 3/2/1996

    Args:
        gasid (int): The gas ID defined in :py:class:`~pyrtlib.atmp.AtmosphericProfiles`
    
    Returns:
        float: The mass of the HITRAN gas ID

    .. note::
        Results are not accurate because amu values need more significant figures.

    """

    if gasid == 0: amus = np.dot(2, 1) + 16
    if gasid == 1: amus = 12 + np.dot(2, 16)
    if gasid == 2: amus = np.dot(3, 16)
    if gasid == 3: amus = np.dot(2, 14) + 16
    if gasid == 4: amus = 12 + 16
    if gasid == 5: amus = 12 + np.dot(4, 1)
    if gasid == 6: amus = np.dot(2, 16)
    if gasid == 7: amus = 14 + 16
    if gasid == 8: amus = 32 + np.dot(2, 16)
    if gasid == 9: amus = 14 + np.dot(2, 16)
    if gasid == 10: amus = 14 + np.dot(3, 1)
    if gasid == 11: amus = 1 + 14 + np.dot(3, 16)
    if gasid == 12: amus = 16 + 1
    if gasid == 13: amus = 1 + 19
    if gasid == 14: amus = 1 + 35
    if gasid == 15: amus = 1 + 80
    if gasid == 16: amus = 1 + 127
    if gasid == 17: amus = 35 + 16
    if gasid == 18: amus = 16 + 12 + 32
    if gasid == 19: amus = np.dot(2, 1) + 12 + 16
    if gasid == 20: amus = 1 + 16 + 35
    if gasid == 21: amus = np.dot(2, 14)
    if gasid == 22: amus = 1 + 12 + 14
    if gasid == 23: amus = 12 + np.dot(3, 1) + 35
    if gasid == 24: amus = np.dot(2, 1) + np.dot(2, 16)
    if gasid == 25: amus = np.dot(2, 12) + np.dot(2, 1)
    if gasid == 26: amus = np.dot(2, 12) + np.dot(6, 1)
    if gasid == 27: amus = 31 + np.dot(3, 1)
    if gasid == 28: amus = 12 + 16 + np.dot(2, 19)
    if gasid == 29: amus = 32 + np.dot(6, 19)
    if gasid == 30: amus = np.dot(2, 1) + 32
    if gasid == 31: amus = 1 + 12 + np.dot(2, 16) + 1
    if gasid == 99: amus = 28.9402753669

    mass_proton = 1.6726485e-27
    mass_molecule = np.dot(mass_proton, amus)

    return mass_molecule


def ppmv2gkg(ppmv: np.ndarray, gasid: int) -> np.ndarray:
    """Convert volume mixing ratio in ppmv to mass mixing ratio in g/kg.

    Args:
        ppmv (np.ndarray): mass mixing ratio (g/kg).
        gasid (int): HITRAN gas id.
    
    Returns:
        np.ndarray: Mass mixing ratio in g/kg

    See also:
        
        :py:meth:`gas_mass` 
    """

    # convert to parts per volume
    ppv = ppmv / 1000000.0
    # multiply by ratio of masses to get mass mixing ratio in g/g
    gg = np.dot(ppv, gas_mass(gasid)) / gas_mass(99)
    # multiply by 1000 to get g/kg
    gkg = np.dot(gg, 1000)

    return gkg


def mr2rh(p: np.ndarray,
          t: np.ndarray,
          w: np.ndarray,
          Tconvert: np.ndarray = None) -> np.ndarray:
    """Determine relative humidity (#) given
    reference pressure (mbar), temperature (t,K), and
    water vapor mass mixing ratio (w,g/kg)
    
    Two RHs are returned: rh1 is with RH defined as the ratio 
    of water vapor partial pressure to saturation vapor pressure and
    rh2 is with RH defined as the ratio of water vapor mixing ratio to 
    saturation mixing ratio.
    
    if input, Tconvert is used as the temperature point to switch
    from using saturation vapor pressure over water to over ice.
    
    DCT 3/5/00


    Args:
        p (np.ndarray): [description].
        t (np.ndarray): [description].
        w (np.ndarray): [description].
        Tconvert (np.ndarray, optional): [description]. Defaults to None.

    Returns:
        np.ndarray: Relative humidity
    """
    # saturation pressure
    esat = satvap(t)
    wsat = satmix(p, t)
    if Tconvert:
        esat = satvap(t, Tconvert)
        wsat = satmix(p, t, Tconvert)

    # H2O partial pressure
    e = mr2e(p, w)
    # RH using ratios of gas pressures
    rh1 = np.dot(100.0, e) / esat
    # RH using WMO definition of relative humidity
    rh2 = np.dot(100, w) / wsat

    return rh1, rh2


def mr2rho(mr: np.ndarray, tk: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Determine water vapor density (g/m3) given
    reference pressure (mbar), temperature (t,K), and
    water vapor mass mixing ratio (g/kg)

    Equations were provided by Holger Linne' from Max Planck Institute.
    
    Nico 2002/05/09 (Looking at rho2mr.email.m from DCT)
    Nico 2018/06/20

    Args:
        mr ([type], optional): [description]. Defaults to None.
        tk ([type], optional): [description]. Defaults to None.
        p ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    rho = np.multiply(np.multiply(mr, p), 0.3477) / (tk)

    # I think the above is an approximation valid within ~1#.
    # To be consistent with Vapor_xxx.m and mr2rh.m (see
    # Compute_Transmittances_for_RTTOV_dsb.m), it should be:
    rvap = np.dot(constants('Rwatvap'), 1e-05)

    eps = 0.621970585
    rho = np.multiply(np.multiply(mr, p), 1.0) / (np.dot(
        (np.dot(1000.0, eps) + mr), rvap)) / (tk)

    return rho


def mr2e(p: np.ndarray, mr: np.ndarray) -> np.ndarray:
    """Compute H2O partial pressure (e,mbar) given
    pressure (p,mbar) and H2O mass mixing ratio (mr,g/kg)

    DCT 3/6/00

    Args:
        p ([type], optional): [description]. Defaults to None.
        mr ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """

    # ratio of water mass to dry air mass
    eps = 0.621970585
    e = np.multiply(p, mr) / (np.dot(1000.0, eps) + mr)

    return e


def e2mr(p: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Compute H2O mass mixing ratio (mr,g/kg) given
    pressure (p,mbar) and H2O partial pressure (e,mbar)
        
    DCT 3/6/00

    Args:
        p ([type], optional): [description]. Defaults to None.
        e ([type], optional): [description]. Defaults to None.
    
    Returns:
        [type]: [description]
    """

    # ratio of water mass to dry air mass
    eps = 0.621970585
    mr = np.multiply(np.dot(1000.0, eps), e) / (p - e)

    return mr


def satmix(p: np.ndarray,
           T: np.ndarray,
           Tconvert: Optional[np.ndarray] = None) -> np.ndarray:
    """Compute saturation mixing ratio [g/kg] given reference pressure, 
    p [mbar] and temperature, T [K].  If Tconvert input, the calculation uses 
    the saturation vapor pressure over ice (opposed to over water) 
    for temperatures less than Tconvert [K].
        
    DCT, updated 3/5/00

    Args:
        p ([type], optional): [description]. Defaults to None.
        T ([type], optional): [description]. Defaults to None.
        Tconvert ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    # warning('off')

    # saturation pressure
    esat = satvap(T)
    if Tconvert:
        esat = satvap(T, Tconvert)

    # saturation mixing ratio
    wsat = e2mr(p, esat)

    return wsat


def satvap(T, Tconvert: Optional[np.ndarray] = None) -> np.ndarray:
    """compute saturation vapor pressure [mbar] given temperature, T [K].
    If Tconvert is input, the calculation uses the saturation vapor 
    pressure over ice (opposed to over water) for temperatures less than 
    Tconvert [K].
        
    DCT, updated 3/5/00

    Args:
        T ([type], optional): [description]. Defaults to None.
        Tconvert ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]

    """
    # saturation pressure over water
    # Goff Gratch formulation, over water
    esat = eswat_goffgratch(T)

    # saturation pressure over ice if needed
    if Tconvert:
        ind = np.nonzero(T <= Tconvert)
        # Goff Gratch formulation, over ice
        esat[ind] = esice_goffgratch(T(ind))

    return esat


def eswat_goffgratch(T: np.ndarray) -> np.ndarray:
    """Compute water vapor saturation pressure over water
    using Goff-Gratch formulation.  Adopted from PvD's 
    svp_water.pro.

    DCT 8/22/00

    Args:
        T ([type], optional): temperature [Kelvin]. Defaults to None.

    Returns:
        svp [type]: saturation pressure [mbar]

    References
    ----------
    .. [1] Goff-Gratch formulation from sixth revised 
            edition of Smithsonian Meteorology Tables.

    .. note::
        svp returned for all values of input T,
        but results not valid for T >= 370 K and 
        T <= 160 K.
    """

    t_sat = 373.16
    t_ratio = t_sat / T
    rt_ratio = 1.0 / t_ratio
    sl_pressure = 1013.246
    c1 = 7.90298
    c2 = 5.02808
    c3 = 1.3816e-07
    c4 = 11.344
    c5 = 0.0081328
    c6 = 3.49149

    tmp = (-1.0 * c1 * (t_ratio - 1.0)) + \
          (c2 * np.log10(t_ratio)) - \
          (c3 * (10.0 ** (c4 * (1.0 - rt_ratio)) - 1.0)) + \
          (c5 * (10.0 ** (-1.0 * c6 * (t_ratio - 1.0)) - 1.0)) + \
          np.log10(sl_pressure)

    svp = 10.0 ** tmp

    return svp


def esice_goffgratch(T: np.ndarray) -> np.ndarray:
    """Compute water vapor saturation pressure over ice
    using Goff-Gratch formulation.  Adopted from PvD's 
    svp_ice.pro.

    DCT 8/22/00

    Args:
        T ([type], optional): [description]. Defaults to None.

    Returns:
        svp [type]: [description]

    References
    ----------
    .. [1] Goff-Gratch formulation from sixth revised 
            edition of Smithsonian Meteorology Tables.

    .. note::
        svp returned for all values of input T,
        but results not valid for T >= 370 K and 
        T <= 160 K.
    """

    ewi = 6.1071
    c1 = 9.09718
    c2 = 3.56654
    c3 = 0.876793
    ratio = 273.15 / T
    tmp = (np.dot(-c1,
                  (ratio - 1.0))) - (np.dot(c2, np.log10(ratio))) + (np.dot(
        c3, (1.0 - (1.0 / ratio)))) + np.log10(ewi)

    svp = 10.0 ** tmp

    return svp


def tk2b_mod(hvk: np.ndarray, T: np.ndarray) -> np.ndarray:
    r"""[summary]

    .. math::

        Btilde=\frac{1}{e^\frac{hvk}{T}-1}

    Args:
        hvk ([type], optional): [description]. Defaults to None.
        T ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]

    .. warning:: add docstring to function
    """
    Btilde = 1.0 / (np.exp(hvk / T) - 1.0)

    return Btilde


def dilec12(f: np.ndarray, tk: np.ndarray) -> np.ndarray:
    """Computes the complex dielectric constant for liquid water,
    with a negative imaginary part representing dissipation.

    Complex logarithm is used here. It should be defined with
    imaginary part in the range -pi to +pi.

    Copyright ? P.W. Rosenkranz  Apr. 15, 2014
    Creative Commons license CC BY-SA

    Args:
        f ([type], optional): frequency in GHz. Defaults to None.
        tk ([type], optional): Kelvin temperature. Defaults to None.

    Returns:
        kappa [type]: complex dielectric constant

    References
    ----------
        .. [1] Static dielectric constant model from Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
        .. [2] Debye term from  W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
        .. [3] B band from P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. v.53(3) pp.1387-93 (2015).


    .. note:: validated for 20<f<220 GHz at 248<tk<273; 1<f<1000 GHz at 273<tk<330.
    """

    tc = tk - 273.15
    z = np.complex(0.0, f)
    theta = 300.0 / tk
    # static dielectric constant model from
    # Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
    kappa = -43.7527 * theta ** 0.05 + 299.504 * theta ** 1.47 - 399.364 * theta ** 2.11 + 221.327 * theta ** 2.31
    # Debye term from
    # W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
    delta = 80.69715 * np.exp(-tc / 226.45)
    sd = 1164.023 * np.exp(-651.4728 / (tc + 133.07))
    kappa = kappa - (delta * z) / (sd + z)
    # B band from
    # P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. v.53(3) pp.1387-93 (2015).
    delta = 4.008724 * np.exp(-tc / 103.05)
    hdelta = delta / 2.0
    f1 = 10.46012 + (0.1454962 * tc) + (0.063267156 * tc ** 2) + (0.00093786645 *
                                                                  tc ** 3)
    # z1 = (-.75,1.) * f1;
    # z2 = (-4500.,2000.)
    z1 = np.complex(-0.75, 1.0) * f1
    z2 = np.complex(-4500.0, 2000.0)
    cnorm = np.log(z2 / z1)
    chip = (hdelta * np.log((z - z2) / (z - z1))) / cnorm
    chij = (hdelta * np.log(
        (z - np.conj(z2)) / (z - np.conj(z1)))) / np.conj(cnorm)
    dchi = chip + chij - delta
    kappa = kappa + dchi

    return kappa


def dcerror(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Csixth-Order Approx To The Complex Error Function Of z=x+iy.
    cerror = exp(-z^2)erfc(-iz)
    This version is double precision and valid in all quadrants.

    P. Rosenkranz  12/11/2018
    2018/12/19 - Nico: first created from dcerror.f

    Args:
        x ([type], optional): [description]. Defaults to None.
        y ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]

    References
    ----------
    .. [1] Hui, Armstrong And Wray, Jqsrt V.19, P.509-516 (1978).
    """

    # IMPLICIT NONE
    # DOUBLE PRECISION X,Y,a(0:6),b(0:6)
    # DOUBLE COMPLEX ASUM,BSUM,ZH,w

    a = np.asarray([
        122.607931777104326, 214.382388694706425, 181.928533092181549,
        93.155580458138441, 30.180142196210589, 5.912626209773153,
        0.564189583562615
    ])

    b = np.asarray([
        122.607931773875350, 352.730625110963558, 457.334478783897737,
        348.703917719495792, 170.354001821091472, 53.992906912940207,
        10.479857114260399
    ])

    # compute w in quadrants 1 or 2
    # from eqs.(13), w(z) = [w(-z*)]*
    # expansion in terms of ZH results in conjugation of w when X changes sign.
    zh = complex(np.abs(y), -x)
    asum = (((((a[6] * zh + a[5]) * zh + a[4]) * zh + a[3]) * zh + a[2]) * zh +
            a[1]) * zh + a[0]
    bsum = (((((
                       (zh + b[6]) * zh + b[5]) * zh + b[4]) * zh + b[3]) * zh + b[2]) * zh +
            b[1]) * zh + b[0]
    w = asum / bsum
    if y >= 0:
        dcerror = w
    else:
        # from eqs.(13), w(z) = 2exp(-z^2)-[w(z*)]*
        dcerror = 2.0 * np.exp(-np.complex(x, y) ** 2) - w

    return dcerror


def pressure_to_height(pressure: float) -> float:
    r"""Convert pressure data to height using the U.S. standard atmosphere [NOAA1976]_.
    The implementation uses the formula outlined in [Hobbs1977]_ pg.60-61.

    .. math:: Z = \frac{T_0}{\Gamma}[1-\frac{p}{p_0}^\frac{R\Gamma}{g}]

    Args:
        pressure (float): The pressure value in mbar

    Returns:
        float: The height in meters to the provided pressure
    """
    t0 = 288.0
    p0 = 1013.25
    gamma = -0.0065
    g = 9.80665
    R = 8.314462618
    Md = 28.96546e-3
    Rd = R / Md

    return (t0 / gamma) * (1 - (pressure / p0) ** (Rd * gamma / g))


def height_to_pressure(height: float) -> float:
    r"""Convert height data to pressures using the U.S. standard atmosphere [NOAA1976]_.
    The implementation inverts the formula outlined in [Hobbs1977]_ pg.60-61.

    .. math:: p = p_0 e^{\frac{g}{R \Gamma} \text{ln}(1-\frac{Z \Gamma}{T_0})}

    Args:
        height (float): The height value in meters

    Returns:
        float: The pressure to the provided height
    """
    t0 = 288.0
    p0 = 1013.25
    gamma = -0.0065
    g = 9.80665
    R = 8.314462618
    Md = 28.96546e-3
    Rd = R / Md

    return p0 * np.exp(
        (g / (Rd * gamma)) * np.log(1 - ((height * gamma) / t0)))


def dewpoint2rh(td: float,
                t: float,
                ice: Optional[bool] = False,
                method: Optional[str] = 'arm') -> float:
    r"""Calculate relative humidity from temperature and dewpoint.
    Value is calculated using the August-Roche-Magnus approximation. [AUGUST]_ [MAGNUS]_.

    .. math:: RH = \frac {\exp(\frac{a T_d}{b+T_d})} {\exp(\frac{a T}{b+T})}

    .. math:: where \ a = 17.625, b = 243.04

    .. math:: RH = \frac {6.1078\times10^{\frac{aT_d}{b + T_db}}}{6.1078\times10^{\frac{a T}{b + T}}}

    .. math:: where \ a = 7.5, b = 265.5

    Args:
        td (float): The dew point temperature in Celsius
        t (float): The temperature in Celsius

    Returns:
        float: The relative humidity to the provided dew point temperature

    References
    ----------
    .. [1] Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus' form approximation of saturation vapor pressure. J. Appl. Meteor., 35, 601–609.
    .. [2] August, E. F., 1828: Ueber die Berechnung der Expansivkraft des Wasserdunstes. Ann. Phys. Chem., 13, 122–137.
    .. [3] Magnus, G., 1844: Versuche über die Spannkräfte des Wasserdampfs. Ann. Phys. Chem., 61, 225–247.
    """
    if method == 'arm':
        a = 17.625
        b = 243.04
        rh = (np.exp((a * td) / (b + td)) / np.exp((a * t) / (b + t)))
    else:
        a = 7.5
        b = 237.3
        if ice:
            a = 9.5
            b = 265.5
        # Vapor pressure (ED):
        ed = 6.1078 * 10 ** ((td * a) / (td + b))

        # saturated vapor pressure (es):
        es = 6.1078 * 10 ** ((t * a) / (t + b))

        # relative humidity = 100 * ed/es
        rh = ed / es

    return rh


def kgkg_to_kgm3(q: np.ndarray, p: np.ndarray, t: np.ndarray) -> np.ndarray:
    r"""Utils function to convert from Kg Kg-1 to kg m-3. 

    NWP models provide cloud liquid and ice water content in units kq kq-1. To convert
    to g m-3 multiply the result of this function to the value in kg kg-1.

    .. math:: LWC = q_{liq} \frac{10^2 P}{R_{moist} T}

    .. math:: R_{moist} = R_{dry} (1 + \frac{1 - \epsilon}{\epsilon} q_{h2o})

    .. math:: \epsilon = \frac{M_{h2o}}{M_{dry}}

    .. math:: R_{dry} = \frac{R}{M_{dry}}

    where:
    :math:`q_{liq}` is the mass mixing ratio for liquid cloud, :math:`P` is the atmospheric pressure 
    in hPa, :math:`T` is the atmospheric temperature in K, :math:`R_{moist}` is the moist 
    air gas constant (in J kg-1 K-1), :math:`R_{dry}` is the gas constant for dry air and :math:`q_{h2o}` is the specific humidity 
    (given as the ratio between the mass of water vapor and the mass of moist air)

    The same equations are used for ice clouds, by replacing LWC by IWC and :math:`q_{liq}` by :math:`q_{ice}`

    Args:
        q (np.ndarray): specific humidity (Kg Kg-1)
        p (np.ndarray): pressure (hPa)
        t (np.ndarray): temperature (K)

    Returns:
        np.ndarray: [description]

    References:
        .. [1] M. Z. Jacobson. Fundamentals of atmospheric modelling. Cambridge Eds., 2005.
    """

    eps = 0.621970585
    rmoist = (constants('Rdry')[0]) * (1 + ((1 - eps) / eps) * q)
    kgm3 = (1e2 * p) / (rmoist * t)

    return kgm3


def ppmv_to_moleculesm3(mr: np.ndarray, p: np.ndarray,
                        t: np.ndarray) -> np.ndarray:
    """For any gas, this function converts mixing ratio  (in ppmv) to number density (molecules/m3).

    Args:
        mr (np.ndarray): mixing ratio in ppmv
        p (np.ndarray): pressure in Pa
        t (np.ndarray): temperature in K

    Returns:
        np.ndarray: [description]
    """
    av = constants('avogadro')[0]
    rg = constants('R')[0]
    n_air = (av * p) / (rg * t)  # (molecules/m3)
    nr_molm3 = n_air * mr * 1e-6

    return nr_molm3


def get_frequencies(instr: Optional[str] = 'hat'):
    """[summary]

    Args:
        instr (str, optional): [description]. Defaults to 'hat'.

    Returns:
        [type]: [description]
    """
    frequencies = {
        'hat': [
            22.2400, 23.0400, 23.8400, 25.4400, 26.2400, 27.8400, 31.4000,
            51.2600, 52.2800, 53.8600, 54.9400, 56.6600, 57.3000, 58.0000
        ],
        'mp3': [
            22.234, 22.500, 23.034, 23.834, 25.000, 26.234, 28.000, 30.000,
            51.248, 51.760, 52.280, 52.804, 53.336, 53.848, 54.400, 54.940,
            55.500, 56.020, 56.660, 57.288, 57.964, 58.800
        ],
        'tem': [
            51.25, 51.75, 52.25, 52.85, 53.35, 53.85, 54.40, 54.90, 55.40,
            56.00, 56.50, 57.00
        ],
        'k2w': [23.8400, 31.4000, 72.5000, 82.5000, 90.0000, 150.000]
    }

    return frequencies.get(instr)

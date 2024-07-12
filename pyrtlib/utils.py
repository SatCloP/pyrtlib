"""
This module contains the utils functions.
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

import types
from typing import Tuple, Optional, Union, List
import sys
from importlib import reload
import numpy as np


def import_lineshape(name: str) -> types.ModuleType:
    """Import a named object from a module in the context of this function.
    Used to import line list for absorption models.

    :meta private:

    Args:
        name (str): Absorption model name.

    Returns:
        types.ModuleType: Dictionary of line list of the absorption model chose.

    See also:
        :py:func:`~pyrtlib.absorption_model.H2OAbsModel.set_ll`

    Examples:
        >>> from pyrtlib.absorption_model import H2OAbsModel
        >>> from pyrtlib.utils import import_lineshape
        >>> H2OAbsModel.model = 'R21SD'
        >>> H2OAbsModel.h2oll = import_lineshape('h2oll')
        >>> H2OAbsModel.h2oll.aself
        array([0., 12.6,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.])
    """
    try:
        module = __import__('pyrtlib._lineshape', globals(),
                            locals(), [name.lower()])
    except ImportError as e:
        return None
    if vars(module)[name.lower()] in sys.modules.values():
        reload(vars(module)[name.lower()])

    return vars(module)[name.lower()]


def constants(name: Optional[str] = None) -> Union[Tuple[float, str], List]:
    """
    This routine will provide values and units for all the
    universal constants necessary to pyrtlib.

    Args:
        name (Optional[str], optional): String specifying which constant is needed. Defaults to None.

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
        ValueError: Raises error wheter no costant available.

    Returns:
        Union[Tuple[float, str], List] Numerical Value of the asked constant and string specifying which units are used

    References
    ----------
    .. [1] P.J. Mohr, B.N. Taylor, and D.B. Newell (2015), "The 2014 CODATA Recommended
            Values of the Fundamental Physical Constants" (Web Version 7.0), http://physics.nist.gov/cuu/index.html
            Values as of 11/12/2015
    .. [2] Janssen, Atmospheric Remote Sensing by Microwave Radiometry, pag.12

    .. note::
        # +/- 0.017 [K] Cosmic Background Temperature (see ref. 2)
        Tcosmicbkg = 2.736;
    """
    constants_dict = {
        'avogadro': [6.022140857e+23, '[mol-1]'],
        'boltzmann': [1.3806579999999998e-23, '[J K-1]'],
        'EarthRadius': [6370.949, '[km]'],
        'gravity': [9.80665, '[m s-2]'],
        'light': [299792458.0, '[m s-1]'],
        'Np2dB': [4.342944819032518, '[dB/Np]'],
        'planck': [6.626075499999999e-34, '[J Hz-1]'],
        'Rdry': [287.04, '[J kg-1 K-1]'],
        'Rwatvap': [461.52, '[J kg-1 K-1]'],
        'Tcosmicbkg': [2.728, '[K]'],
        'R': [8.31446261815324, '[J mol-1 K-1]']
    }

    if not name:
        cs = list(constants_dict.keys())
    else:
        try:
            cs = constants_dict[name]
        except KeyError as e:
            raise ValueError(
                'No constant avalaible with this name: {} . Type constants() without argument to get the list of constants'.format(
                    name))

    return cs


def gas_mass(gasid: int) -> float:
    """
    Returns the mass of the HITRAN gas ID

    Args:
        gasid (int): The gas ID defined in :py:class:`~pyrtlib.climatology.AtmosphericProfiles`

    Returns:
        float: The mass of the HITRAN gas ID.
    """

    if gasid == 0:
        amus = np.dot(2, 1) + 16
    if gasid == 1:
        amus = 12 + np.dot(2, 16)
    if gasid == 2:
        amus = np.dot(3, 16)
    if gasid == 3:
        amus = np.dot(2, 14) + 16
    if gasid == 4:
        amus = 12 + 16
    if gasid == 5:
        amus = 12 + np.dot(4, 1)
    if gasid == 6:
        amus = np.dot(2, 16)
    if gasid == 7:
        amus = 14 + 16
    if gasid == 8:
        amus = 32 + np.dot(2, 16)
    if gasid == 9:
        amus = 14 + np.dot(2, 16)
    if gasid == 10:
        amus = 14 + np.dot(3, 1)
    if gasid == 11:
        amus = 1 + 14 + np.dot(3, 16)
    if gasid == 12:
        amus = 16 + 1
    if gasid == 13:
        amus = 1 + 19
    if gasid == 14:
        amus = 1 + 35
    if gasid == 15:
        amus = 1 + 80
    if gasid == 16:
        amus = 1 + 127
    if gasid == 17:
        amus = 35 + 16
    if gasid == 18:
        amus = 16 + 12 + 32
    if gasid == 19:
        amus = np.dot(2, 1) + 12 + 16
    if gasid == 20:
        amus = 1 + 16 + 35
    if gasid == 21:
        amus = np.dot(2, 14)
    if gasid == 22:
        amus = 1 + 12 + 14
    if gasid == 23:
        amus = 12 + np.dot(3, 1) + 35
    if gasid == 24:
        amus = np.dot(2, 1) + np.dot(2, 16)
    if gasid == 25:
        amus = np.dot(2, 12) + np.dot(2, 1)
    if gasid == 26:
        amus = np.dot(2, 12) + np.dot(6, 1)
    if gasid == 27:
        amus = 31 + np.dot(3, 1)
    if gasid == 28:
        amus = 12 + 16 + np.dot(2, 19)
    if gasid == 29:
        amus = 32 + np.dot(6, 19)
    if gasid == 30:
        amus = np.dot(2, 1) + 32
    if gasid == 31:
        amus = 1 + 12 + np.dot(2, 16) + 1
    if gasid == 99:
        amus = 28.9402753669

    mass_proton = 1.6726485e-27
    mass_molecule = np.dot(mass_proton, amus)

    return mass_molecule


def ppmv2gkg(ppmv: np.ndarray, gasid: int) -> np.ndarray:
    """Convert volume mixing ratio in ppmv to mass mixing ratio in g/kg.

    Args:
        ppmv (numpy.ndarray): mass mixing ratio (g/kg).
        gasid (int): HITRAN gas id.

    Returns:
        numpy.ndarray: Mass mixing ratio in g/kg

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
          Tconvert: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """Determine relative humidity (rh) given reference pressure (mbar), temperature (K), and
    water vapor mass mixing ratio (g/kg)

    Two RHs are returned: rh1 is with RH defined as the ratio of water vapor partial pressure
    to saturation vapor pressure and rh2 is with RH defined as the ratio of water vapor mixing ratio to
    saturation mixing ratio.

    if input, Tconvert is used as the temperature point to switch
    from using saturation vapor pressure over water to over ice.

    Args:
        p (numpy.ndarray): Pressure profile (mb).
        t (numpy.ndarray): Temperature profile (K).
        w (numpy.ndarray): Water Vapor Mixing ratio (g/kg).
        Tconvert (numpy.ndarray, optional): Threshold temperature below which saturation 
            water pressure is calculated over ice instead of liquid water. Defaults to None.

    Returns:
        numpy.ndarray: Relative humidity using ratios of gas pressures (%)
        numpy.ndarray: Relative humidity using WMO definition (%)
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


def mr2rho(mr: np.ndarray, t: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Determine water vapor density (:math:`g/m^3`) given reference pressure (mbar), temperature (K), and
    water vapor mass mixing ratio (g/kg)

    Equations were provided by Holger Linne' from Max Planck Institute.

    Args:
        mr (numpy.ndarray): Mixing ratio (g/kg)
        t (numpy.ndarray): Temperature profiles (K)
        p (numpy.ndarray): Pressure profiles (mb).

    Returns:
        numpy.ndarray: Water Vapor Density (:math:`g/m^3`)
    """

    rho = np.multiply(np.multiply(mr, p), 0.3477) / t

    # I think the above is an approximation valid within ~1#.
    # To be consistent with Vapor_xxx.m and mr2rh.m (see
    # Compute_Transmittances_for_RTTOV_dsb.m), it should be:
    rvap = np.dot(constants('Rwatvap')[0], 1e-05)

    eps = 0.621970585
    rho = np.multiply(np.multiply(mr, p), 1.0) / (np.dot(
        (np.dot(1000.0, eps) + mr), rvap)) / t

    return rho


def mr2e(p: np.ndarray, mr: np.ndarray) -> np.ndarray:
    """Compute :math:`H_2O` partial pressure (mbar) given pressure (mbar)
    and :math:`H_2O` mass mixing ratio (g/kg)

    DCT 3/6/00

    Args:
        p (numpy.ndarray): Pressure profile (mb).
        mr (numpy.ndarray): Mixing Ratio (g/kg).

    Returns:
        numpy.ndarray: :math:`H_2O` partial pressure (mb).
    """

    # ratio of water mass to dry air mass
    eps = 0.621970585
    e = np.multiply(p, mr) / (np.dot(1000.0, eps) + mr)

    return e


def rho2rh(rho: np.ndarray, t: np.ndarray, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Convert water vapor density to relative humidity.

    Args:
        rho (np.ndarray): Water vapor density (:math:`g/m^3`)
        t (np.ndarray): Temperature (K)
        p (np.ndarray): Pressure (mb).

    Returns:
        numpy.ndarray: Relative humidity using ratios of gas pressures (%)
        numpy.ndarray: Relative humidity using WMO definition (%)
    """
    e = (rho * t)/216.7
    mr = e2mr(p, e)

    rh1, rh2 = mr2rh(p, t, mr)

    return rh1, rh2


def rho2mr(rho: np.ndarray, t: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Determine water vapor mass mixing ratio (g/kg) given reference pressure (mbar), 
    temperature (t,K), and water vapor density (:math:`g/m^3`).

    Args:
        rho (np.ndarray): Water vapor density (:math:`g/m^3`)
        t (np.ndarray): Temperature (K)
        p (np.ndarray): Pressure (mb).

    Returns:
        np.ndarray: H2O Mass Mixing Ratio (g/kg)
    """

    mr = rho * t / (p * 0.3477)

    rvap = constants('Rwatvap')[0] * 1e-5  # [J kg-1 K-1] -> [hPa * m2 g-1 K-1]
    eps = 0.621970585  # Rdry/Rvap
    mr = rho * (1e3*eps*rvap) / (p/t - rvap * rho)

    return mr


def e2mr(p: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Compute :math:`H_2O` mass mixing ratio (g/kg) given pressure (mbar)
    and :math:`H_2O` partial pressure (mbar)

    Args:
        p (numpy.ndarray): Pressure (mb).
        e (numpy.ndarray): :math:`H_2O` partial pressure (mb).

    Returns:
        numpy.ndarray: :math:`H_2O` Mass Mixing Ratio (g/kg)
    """

    # ratio of water mass to dry air mass
    eps = 0.621970585
    mr = np.multiply(np.dot(1000.0, eps), e) / (p - e)

    return mr


def satmix(p: np.ndarray,
           t: np.ndarray,
           Tconvert: Optional[np.ndarray] = None) -> np.ndarray:
    """Compute saturation mixing ratio (g/kg) given reference pressure,
    p (mbar]) and temperature, T (K).  If Tconvert input, the calculation uses
    the saturation vapor pressure over ice (opposed to over water)
    for temperatures less than Tconvert [K].

    DCT, updated 3/5/00

    Args:
        p (numpy.ndarray): Pressure profile (mb).
        t (numpy.ndarray): Temperature profile (K).
        Tconvert (Optional[numpy.ndarray], optional): Threshold temperature below which saturation 
            water pressure is calculated over ice instead of liquid water. Defaults to None.

    Returns:
        numpy.ndarray: Saturation mixing ratio (g/kg).
    """
    # warning('off')

    # saturation pressure
    esat = satvap(t)
    if Tconvert:
        esat = satvap(t, Tconvert)

    # saturation mixing ratio
    wsat = e2mr(p, esat)

    return wsat


def satvap(t: np.ndarray, Tconvert: Optional[np.ndarray] = None) -> np.ndarray:
    """Compute saturation vapor pressure (mbar) given temperature, T (K).
    If Tconvert is input, the calculation uses the saturation vapor
    pressure over ice (opposed to over water) for temperatures less than
    Tconvert (K).

    Args:
        t (numpy.ndarray): Temperature profile (K).
        Tconvert (Optional[numpy.ndarray], optional): Threshold temperature below which saturation 
            water pressure is calculated over ice instead of liquid water. Defaults to None.

    Returns:
        numpy.ndarray: Saturation vapor pressure (mbar).

    """
    # saturation pressure over water
    # Goff Gratch formulation, over water
    esat = eswat_goffgratch(t)

    # saturation pressure over ice if needed
    if Tconvert:
        ind = np.nonzero(t <= Tconvert)
        # Goff Gratch formulation, over ice
        esat[ind] = esice_goffgratch(t[ind])

    return esat


def eswat_goffgratch(t: np.ndarray) -> np.ndarray:
    """Compute water vapor saturation pressure over water using Goff-Gratch formulation.

    Args:
        t (numpy.ndarray): Temperature profile (K).

    Returns:
        numpy.ndarray: Water vapor saturation pressure over water (mb).

    References
    ----------
    .. [1] Goff-Gratch formulation from sixth revised edition of Smithsonian Meteorology Tables.

    .. note::
        svp returned for all values of input T, but results not valid for T >= 370 K and T <= 160 K.
    """

    t_sat = 373.16
    t_ratio = t_sat / t
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


def esice_goffgratch(t: np.ndarray) -> np.ndarray:
    """Compute water vapor saturation pressure over ice using Goff-Gratch formulation.

    Args:
        t (numpy.ndarray): Temperature profile (K).

    Returns:
        numpy.ndarray: Water vapor saturation pressure over ice (mb).

    References
    ----------
    .. [1] Goff-Gratch formulation from sixth revised edition of Smithsonian Meteorology Tables.

    .. note::
        svp returned for all values of input T, but results not valid for T >= 370 K and T <= 160 K.
    """

    ewi = 6.1071
    c1 = 9.09718
    c2 = 3.56654
    c3 = 0.876793
    ratio = 273.15 / t
    tmp = (np.dot(-c1,
                  (ratio - 1.0))) - (np.dot(c2, np.log10(ratio))) + (np.dot(
                      c3, (1.0 - (1.0 / ratio)))) + np.log10(ewi)

    svp = 10.0 ** tmp

    return svp


def tk2b_mod(hvk: np.ndarray, t: np.ndarray) -> np.ndarray:
    r"""Get modified Planck function (Planck function without the constants :math:`\frac{2h\nu^3}{c^2}`)
    by T and hvk (Planck constant * frequency) / Boltzmann constant, (equation (4) from [Schroeder-Westwater-1991]_)

    .. math:: \tilde{B} = \frac{1}{ e^{\frac{h\nu}{k_{B}T}}-1}

    Args:
        hvk (numpy.ndarray): (Planck constant * frequency) / Boltzmann constant.
        t (numpy.ndarray): Temperature (K)

    Returns:
        numpy.ndarray: Modified Planck function.
    """
    Btilde = 1.0 / (np.exp(hvk / t) - 1.0)

    return Btilde


def dilec12(f: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Computes the complex dielectric constant for liquid water, with a
    negative imaginary part representing dissipation.

    Complex logarithm is used here. It should be defined with
    imaginary part in the range -pi to +pi.

    Args:
        f (numpy.ndarray): Frequency (GHz)
        t (numpy.ndarray): Temeprature (K)

    Returns:
        numpy.ndarray: Dielectric constant for liquid water.

    References
    ----------
        .. [1] Static dielectric constant model from Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
        .. [2] Debye term from  W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
        .. [3] B band from P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. v.53(3) pp.1387-93 (2015).


    .. note:: validated for 20<f<220 GHz at 248<tk<273; 1<f<1000 GHz at 273<tk<330.
    """

    tc = t - 273.15
    z = complex(0.0, f)
    theta = 300.0 / t
    # static dielectric constant model from
    # Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
    kappa = -43.7527 * theta ** 0.05 + 299.504 * theta ** 1.47 - \
        399.364 * theta ** 2.11 + 221.327 * theta ** 2.31
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
    z1 = complex(-0.75, 1.0) * f1
    z2 = complex(-4500.0, 2000.0)
    cnorm = np.log(z2 / z1)
    chip = (hdelta * np.log((z - z2) / (z - z1))) / cnorm
    chij = (hdelta * np.log(
        (z - np.conj(z2)) / (z - np.conj(z1)))) / np.conj(cnorm)
    dchi = chip + chij - delta
    kappa = kappa + dchi

    return kappa


def _dcerror(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    r"""Sixth-Order Approx To The Complex Error Function of

    .. math:: z = x+iy
    .. math:: cerror = \exp(z^2)\times erfc(-iz)

    This version is double precision and valid in all quadrants.

    Args:
        x (numpy.ndarray): _description_
        y (numpy.ndarray): _description_

    Returns:
        numpy.ndarray: _description_

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
        dcerror = 2.0 * np.exp(-complex(x, y) ** 2) - w

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


def virtual_temperature(t: np.ndarray, mr: np.ndarray) -> np.ndarray:
    r"""Calculate virtual temperature.
    This calculation must be given an air parcel's temperature and mixing ratio.
    The implementation uses the formula outlined in [Hobbs2006]_ pg.80.

    .. math:: T_v = T \frac{\text{w} + \epsilon}{\epsilon\,(1 + \text{w})}

    Args:
        t (numpy.ndarray): Air temperature (K)
        mr (numpy.ndarray): Mass mixing ratio (dimensionless kg/kg-1)

    Returns:
        numpy.ndarray: Corresponding virtual temperature of the parcel

    Examples:
        >>> from pyrtlib.utils import virtual_temperature
        >>> virtual_temperature(283.2, 12*1e-3)
        285.2412547754703

    .. note::
        This function is based on metpy.calc.virtual_temperature method.
    """

    molecular_weight_ratio = constants("Rdry")[0] / constants("Rwatvap")[0]

    return t * ((mr + molecular_weight_ratio)
                / (molecular_weight_ratio * (1 + mr)))


def _thickness_hydrostatic(p: np.ndarray, t: np.ndarray, mr: Optional[np.ndarray] = None) -> np.float32:
    """
    :meta private:
    """

    R = 8.314462618
    Md = 28.96546e-3
    Rd = R / Md
    g = 9.80665
    if mr is None:
        layer_p, layer_virttemp = p, t
    else:
        layer_p = p
        layer_virttemp = virtual_temperature(t, mr)

    return (
        -Rd / g * np.trapz(layer_virttemp, np.log(layer_p))
    )


def atmospheric_tickness(p: np.ndarray, t: np.ndarray, mr: Optional[np.ndarray] = None) -> np.ndarray:
    r"""Calculate the thickness of a layer via the hypsometric equation.
    This thickness calculation uses the pressure and temperature profiles (and optionally
    mixing ratio) via the hypsometric equation with virtual temperature adjustment.

    .. math:: Z_2 - Z_1 = -\frac{R_d}{g} \int_{p_1}^{p_2} T_v d\ln p,

    Which is based off of Equation 3.24 in [Hobbs2006]_.
    This assumes a hydrostatic atmosphere.

    Args:
        p (numpy.ndarray): Atmospheric pressure profile (mb)
        t (numpy.ndarray): Atmospheric temperature profile (K).
        mr (Optional[numpy.ndarray], optional): Mass mixing ratio (dimensionless kg/kg-1). Defaults to None.

    Returns:
        numpy.ndarray: The thickness of the layers in kilometers

    .. note::
        This function is based on metpy.calc.thickness_hydrostatic method.
    """

    h = np.zeros(p.shape)
    for i in range(p.shape[0]):
        if i == 0:
            continue
        hi = _thickness_hydrostatic(
            p[0:i+1], t[0:i+1], mr[0:i+1] if mr is not None else None)
        h[i] = hi/1e3

    return h


def dewpoint2rh(td: float,
                t: float,
                ice: Optional[bool] = False,
                method: Optional[str] = 'arm') -> float:
    r"""Calculate relative humidity from temperature and dewpoint.
    Value is calculated using the August-Roche-Magnus approximation. [August-1828]_ [Magnus-1844]_.

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
    .. [1] [Alduchov-1996]_
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
    r"""Utils function to convert from :math:`kg/kg` to :math:`kg/m^3`. [Jacobson]_

    NWP models provide cloud liquid and ice water content in units :math:`kg/kg`. To convert
    to :math:`g/m^3` multiply the result of this function to the value in :math:`kg/kg`.

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
        q (numpy.ndarray): specific humidity (:math:`kg/kg`)
        p (numpy.ndarray): pressure (hPa)
        t (numpy.ndarray): temperature (K)

    Returns:
        numpy.ndarray: [description]

    References:
        .. [1] M.Z., Jacobson. Fundamentals of atmospheric modelling. Cambridge Eds., 2005.
    """

    eps = 0.621970585
    rmoist = (constants('Rdry')[0]) * (1 + ((1 - eps) / eps) * q)
    kgm3 = (1e2 * p) / (rmoist * t)

    return kgm3


def ppmv_to_moleculesm3(mr: np.ndarray, p: np.ndarray,
                        t: np.ndarray) -> np.ndarray:
    """For any gas, this function converts mixing ratio  (in ppmv) to number density (:math:`molecules/m^3`).

    Args:
        mr (numpy.ndarray): mixing ratio in ppmv
        p (numpy.ndarray): pressure in Pa
        t (numpy.ndarray): temperature in K

    Returns:
        numpy.ndarray: [description]
    """
    av = constants('avogadro')[0]
    rg = constants('R')[0]
    n_air = (av * p) / (rg * t)  # (molecules/m3)
    nr_molm3 = n_air * mr * 1e-6

    return nr_molm3


def get_frequencies(instrument: Optional[str] = 'hat') -> List:
    """Get frequencies list from main ground instrument.

    Args:
        instrument (Optional[str], optional): Abbreviation of radiometer name. Defaults to 'hat'.

    Returns:
        List: Frequencies list (GHz) of the radiometer chosen.
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
        'arm': [22.235, 23.035, 23.835, 26.235, 30.000, 51.250, 52.280, 53.850, 54.940,
                56.660, 57.290, 58.800],
        'tem': [
            51.25, 51.75, 52.25, 52.85, 53.35, 53.85, 54.40, 54.90, 55.40,
            56.00, 56.50, 57.00
        ],
        'k2w': [23.8400, 31.4000, 72.5000, 82.5000, 90.0000, 150.000]
    }

    try:
        return np.array(frequencies[instrument.lower()])
    except KeyError:
        raise ValueError(
            f"Invalid instrument name. Available instruments are: {list(frequencies.keys())}")


def to_kelvin(t: np.ndarray) -> np.ndarray:
    """Convert T from Celsius to Kelvin

    Args:
        t (numpy.ndarray): Temperature (°C)

    Returns:
        numpy.ndarray: Temperature (K)
    """
    t_k = t + 273.25

    return t_k


def to_celsius(t: np.ndarray) -> np.ndarray:
    """Convert T from Kelvin to Celsius

    Args:
        t (numpy.ndarray): Temperature (K)

    Returns:
        numpy.ndarray: Temperature (°C)
    """
    t_c = t - 273.25

    return t_c


def get_frequencies_sat(instrument: str) -> np.ndarray:
    """Get frequencies list from main satellite sensors.

    Args:
        instrument (str): Instrument from which getting frequencies

    Returns:
        numpy.ndarray: Frequencies (GHz) of the instrument selected

    See Also:
        :py:meth:`pyrtlib.utils.get_frequencies`

    Example:
        .. code-block:: python

            from pyrtlib.utils import get_frequencies_sat
            mwi = get_frequencies_sat("MWI")
    """
    cf = 183.31
    cf2 = 243.20
    cf3 = 325.15
    cf4 = 448.00
    cf6 = 664.00
    cf118 = 118.7503
    cf165 = 165.5
    cf53 = 57.290344
    cf57 = 57.290344

    instrument = instrument.upper()

    if instrument == 'SAPHIR':
        freq = np.array([cf-11, cf-6.8, cf-4.2, cf-2.8, cf-1.1, cf-0.2,
                        cf+0.2, cf+1.1, cf+2.8, cf+4.2, cf+6.8, cf+11])

    elif instrument == 'AMSU':
        freq = np.array([89.0, 150.0, cf-7.0, cf-3.0, cf-1.0,
                        cf+1.0, cf+3.0, cf+7.0])  # AMSU-B
    elif instrument == 'ATMS':
        freq = np.array([88.2, 165.5, cf-7.0, cf-4.5, cf-3.0, cf -
                        1.8, cf-1.0, cf+1.0, cf+1.80, cf+3.0, cf+4.5, cf+7.0])
    elif instrument == 'MWI':
        freq = np.array([18.7, 23.8, 31.4, 50.3, 52.61, 53.24, 53.75, 89,
                        cf118-3.2, cf118-2.1, cf118-1.4, cf118 -
                         1.2, cf118+1.2, cf118+1.4, cf118+2.1, cf118+3.2,
                        cf165-0.725, cf165+0.725,
                        cf-8.4, cf-6.1, cf-4.9, cf-3.4, cf-2.0, cf+2.0, cf+3.4, cf+4.9, cf+6.1, cf+8.4])
    elif instrument == 'MWS':
        freq = np.array([23.8, 31.4, 50.3, 52.8,
                        53.246-0.08, 53.246+0.08, 53.596-0.115, 53.596+0.115, 53.948-0.081, 53.948+0.081,
                        54.4, 54.94, 55.5, cf57,
                        cf57-0.217, cf57-0.217,
                        cf57-0.3222-0.048, cf57-0.3222+0.048, cf57+0.3222-0.048, cf57+0.3222+0.048,
                        cf57-0.3222-0.022, cf57-0.3222+0.022, cf57+0.3222-0.022, cf57+0.3222+0.022,
                        cf57-0.3222-0.010, cf57-0.3222+0.010, cf57+0.3222-0.010, cf57+0.3222+0.010,
                        cf57-0.3222-0.045, cf57-0.3222+0.045, cf57+0.3222-0.045, cf57+0.3222+0.045,
                        89,
                        cf165-0.725, cf165+0.725,
                        cf-7.0, cf-4.5, cf-3.0, cf-1.8, cf-1.0, cf+1.0, cf+1.8, cf+3.0, cf+4.5, cf+7.0,
                        229])
    elif instrument == 'ICI':
        freq = np.array([cf-7.0, cf-3.4, cf-2.0, cf-0.2, cf+0.2, cf+2.0, cf+3.4, cf+7.0,
                        cf2-2.5, cf2+2.5,
                        cf3-9.5, cf3-3.5, cf3-1.5, cf3+1.5, cf3+3.5, cf3+9.5,
                        cf4-7.2, cf4-3.0, cf4-1.4, cf4+1.4, cf4+3.0, cf4+7.2,
                        cf6-4.2, cf6+4.2])
    else:
        raise ValueError("Instrument not found")

    return freq

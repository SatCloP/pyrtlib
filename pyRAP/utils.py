# -*- coding: utf-8 -*-
"""
***************************************************************************
    utils.py
    ---------------------
    Date            : March 2021
    Copyright       : (C) 2021 by CNR-IMAA
    Email           : 
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = ''
__date__ = 'March 2021'
__copyright__ = '(C) 2021, CNR-IMAA'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from copy import copy
import numpy as np


def function(f):
    def helper(*args, **kwargs):
        helper.nargin = len(args)
        # helper.varargin = cellarray(args)
        return f(*args, **kwargs)

    return helper


@function
def constants(string=None, *args, **kwargs):
    """This routine will provide values and units for all the
    universal constants that I needed in my work.

    Nico, 2000

    History
    2003/06/06 - Added 'Rdry' and 'Rwatvap'
    2005/12/10 - Added 'Tcosmicbkg'
    2021/03/09 - Added 'avogadro', 'gravity'

    Args:
        string (str, optional): String specifying which constant is needed. Defaults to None.
        Avalaible right now:
              'avogadro'          Avogadro number [mol-1]
              'boltzmann'         Boltzmann constant [J K-1]
              'EarthRadius'       Earth radius [km]
              'light'             Light speed [m s-1]
              'Np2dB'             Neper to Decibel [dB/Np]
              'planck'            Planck constant [J Hz-1]
              'Rdry'              Gas constant of dry air [J kg-1 K-1]
              'Rwatvap'           Gas constant of water vapor [J kg-1 K-1]
              'Tcosmicbkg'        Cosmic Background Temperature [K]
    Raises:
        ValueError: [description]

    Returns:
        float const: Numerical Value of the asked constant
        str units: String specifying which units are used

    Reference (not for all):
    ----------
    .. [1] P.J. Mohr, B.N. Taylor, and D.B. Newell (2015), "The 2014 CODATA Recommended
            Values of the Fundamental Physical Constants" (Web Version 7.0), http://physics.nist.gov/cuu/index.html
            Values as of 11/12/2015
    """

    if string == 'avogadro':
        nA = 6.022140857e+23
        units = '[mol-1]'
        out = copy(nA)
    elif string == 'boltzmann':
        K = np.dot(1.380658, 1e-23)
        units = '[J K-1]'
        out = copy(K)
    elif string == 'EarthRadius':
        R = 6370.949
        units = '[Km]'
        out = copy(R)
    elif string == 'gravity':
        g = 9.80665
        units = '[m s-2]'
        out = copy(g)
    elif string == 'light':
        c = 299792458
        units = '[m s-1]'
        out = copy(c)
    elif string == 'Np2dB':
        Np2dB = np.dot(10, log10(exp(1)))
        units = '[dB/Np]'
        out = copy(Np2dB)
    elif string == 'planck':
        h = np.dot(6.6260755, 1e-34)
        units = '[J Hz-1]'
        out = copy(h)
    elif string == 'Rdry':
        Rd = 287.04
        units = '[J kg-1 K-1]'
        out = copy(Rd)
    elif string == 'Rwatvap':
        Rv = 461.5
        units = '[J kg-1 K-1]'
        out = copy(Rv)
    elif string == 'Tcosmicbkg':
        # Tcos = 2.736; # +/- 0.017 [K] Cosmic Background Temperature, 
        # from Janssen, Atmospheric Remote Sensing by Microwave Radiometry, pag.12
        Tcos = 2.728
        units = '[K]'
        out = copy(Tcos)
    else:
        raise ValueError('No constant avalaible with this name: {} . Sorry...'.format(strng))

    return out, units


@function
def gas_mass(gasid=None, *args, **kwargs):
    """Returns the mass of the HITRAN gas ID gasid
    DCT 3/2/1996

    Args:
        gasid ([type], optional): [description]. Defaults to None.
    
    Returns:
        [type]: [description]

    .. note::
        Results are not accurate because amu values need more significant figures
    """

    # if gasid ==  1; H20
    # if gasid ==  2; CO2
    # if gasid ==  3; 03
    # if gasid ==  4; N2O
    # if gasid ==  5; CO
    # if gasid ==  6; CH4
    # if gasid ==  7; O2
    # if gasid ==  8; NO
    # if gasid ==  9; SO2
    # if gasid == 10; NO2
    # if gasid == 11; NH3
    # if gasid == 12; HNO3
    # if gasid == 13; OH
    # if gasid == 14; HF
    # if gasid == 15; HCL
    # if gasid == 16; HBR
    # if gasid == 17; HI
    # if gasid == 18; CLO
    # if gasid == 19; OCS
    # if gasid == 20; H2CO
    # if gasid == 21; HOCL
    # if gasid == 22; N2
    # if gasid == 23; HCN
    # if gasid == 24; CH3CL
    # if gasid == 25; H2O2
    # if gasid == 26; C2H2
    # if gasid == 27; C2H6
    # if gasid == 28; PH3
    # if gasid == 29; COF2
    # if gasid == 30; SF6
    # if gasid == 31; H2S
    # if gasid == 32; HCOOH
    # if gasid == 99; AIR

    if gasid == 1: amus = np.dot(2, 1) + 16
    if gasid == 2: amus = 12 + np.dot(2, 16)
    if gasid == 3: amus = np.dot(3, 16)
    if gasid == 4: amus = np.dot(2, 14) + 16
    if gasid == 5: amus = 12 + 16
    if gasid == 6: amus = 12 + np.dot(4, 1)
    if gasid == 7: amus = np.dot(2, 16)
    if gasid == 8: amus = 14 + 16
    if gasid == 9: amus = 32 + np.dot(2, 16)
    if gasid == 10: amus = 14 + np.dot(2, 16)
    if gasid == 11: amus = 14 + np.dot(3, 1)
    if gasid == 12: amus = 1 + 14 + np.dot(3, 16)
    if gasid == 13: amus = 16 + 1
    if gasid == 14: amus = 1 + 19
    if gasid == 15: amus = 1 + 35
    if gasid == 16: amus = 1 + 80
    if gasid == 17: amus = 1 + 127
    if gasid == 18: amus = 35 + 16
    if gasid == 19: amus = 16 + 12 + 32
    if gasid == 20: amus = np.dot(2, 1) + 12 + 16
    if gasid == 21: amus = 1 + 16 + 35
    if gasid == 22: amus = np.dot(2, 14)
    if gasid == 23: amus = 1 + 12 + 14
    if gasid == 24: amus = 12 + np.dot(3, 1) + 35
    if gasid == 25: amus = np.dot(2, 1) + np.dot(2, 16)
    if gasid == 26: amus = np.dot(2, 12) + np.dot(2, 1)
    if gasid == 27: amus = np.dot(2, 12) + np.dot(6, 1)
    if gasid == 28: amus = 31 + np.dot(3, 1)
    if gasid == 29: amus = 12 + 16 + np.dot(2, 19)
    if gasid == 30: amus = 32 + np.dot(6, 19)
    if gasid == 31: amus = np.dot(2, 1) + 32
    if gasid == 32: amus = 1 + 12 + np.dot(2, 16) + 1
    if gasid == 99: amus = 28.9402753669

    mass_proton = 1.6726485e-27
    mass_molecule = np.dot(mass_proton, amus)

    return mass_molecule


@function
def ppmv2gkg(ppmv=None, gasid=None, *args, **kwargs):
    """Convert volume mixing ratio in ppmv to mass mixing ratio in g/kg.

    Args:
        ppmv ([type], optional): mass mixing ratio (g/kg). Defaults to None.
        gasid ([type], optional): HITRAN gas id. Defaults to None.
    
    Returns:
        [type]: [description]

    See also
    --------
    gass_mass
    """

    # convert to parts per volume
    ppv = ppmv / 1000000.0
    # multiply by ratio of masses to get mass mixing ratio in g/g
    gg = np.dot(ppv, gas_mass(gasid)) / gas_mass(99)
    # multiply by 1000 to get g/kg
    gkg = np.dot(gg, 1000)

    return gkg


@function
def mr2rh(p=None, t=None, w=None, Tconvert=None, *args, **kwargs):
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
        p ([type], optional): [description]. Defaults to None.
        t ([type], optional): [description]. Defaults to None.
        w ([type], optional): [description]. Defaults to None.
        Tconvert ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    nargin = mr2rh.nargin

    # saturation pressure
    if nargin == 3:
        esat = satvap(t)
        wsat = satmix(p, t)
    if nargin == 4:
        esat = satvap(t, Tconvert)
        wsat = satmix(p, t, Tconvert)

    # H2O partial pressure
    e = mr2e(p, w)
    # RH using ratios of gas pressures
    rh1 = np.dot(100.0, e) / esat
    # RH using WMO definition of relative humidity
    rh2 = np.dot(100, w) / wsat

    return rh1, rh2


@function
def mr2rho(mr=None, tk=None, p=None, *args, **kwargs):
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


@function
def mr2e(p=None, mr=None, *args, **kwargs):
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


@function
def e2mr(p=None, e=None, *args, **kwargs):
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


@function
def satmix(p=None, T=None, Tconvert=None, *args, **kwargs):
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

    nargin = satmix.nargin

    # warning('off')

    # saturation pressure
    if nargin == 2:
        esat = satvap(T)
    if nargin == 3:
        esat = satvap(T, Tconvert)

    # saturation mixing ratio
    wsat = e2mr(p, esat)

    return wsat


@function
def satvap(T=None, Tconvert=None, *args, **kwargs):
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
    nargin = satvap.nargin

    # saturation pressure over water
    # Goff Gratch formulation, over water
    esat = eswat_goffgratch(T)

    # saturation pressure over ice if needed
    if nargin == 2:
        ind = find(T <= Tconvert)
        # Goff Gratch formulation, over ice
        esat[ind] = esice_goffgratch(T(ind))

    return esat


@function
def eswat_goffgratch(T=None, *args, **kwargs):
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
    tmp = (np.dot(np.dot(-1.0, c1), (t_ratio - 1.0))) + (np.dot(c2, np.log10(t_ratio))) - (
        np.dot(c3, (10.0 ** (np.dot(c4, (1.0 - rt_ratio))) - 1.0))) + (np.dot(
        c5, (10.0 ** (np.dot(np.dot(-1.0, c6),
                             (t_ratio - 1.0))) - 1.0))) + np.log10(sl_pressure)
    svp = 10.0 ** tmp

    return svp


@function
def esice_goffgratch(T=None, *args, **kwargs):
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
    tmp = (np.dot(-c1, (ratio - 1.0))) - (np.dot(c2, np.log10(ratio))) + (np.dot(
        c3, (1.0 - (1.0 / ratio)))) + np.log10(ewi)

    svp = 10.0 ** tmp

    return svp

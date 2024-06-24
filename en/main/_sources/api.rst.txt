==============
API references
==============

Main class
==========

The main class which computes brightness temperatures (Tb), mean radiating temperature (Tmr), and integrated absorption (Tau) for 
clear or cloudy conditions. Also returns all integrated quantities that the original TBMODEL, Cyber Version, returned ([Schroeder-Westwater-1991]_).

.. autosummary::
    :toctree: generated/

    pyrtlib.tb_spectrum.TbCloudRTE


Example:

Compute downwelling (:code:`rte.satellite == False`) brightness temperature for a typical Tropical Atmosphere.

.. plot::
    :include-source: true

    from pyrtlib.tb_spectrum import TbCloudRTE
    from pyrtlib.climatology import AtmosphericProfiles as atmp
    from pyrtlib.utils import ppmv2gkg, mr2rh

    z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100

    ang = np.array([90.])
    frq = np.arange(20, 201, 1)

    rte = TbCloudRTE(z, p, t, rh, frq, ang)
    rte.init_absmdl('R19SD')
    rte.satellite = False
    df = rte.execute()
    df = df.set_index(frq)
    df.tbtotal.plot(figsize=(12,8), xlabel="Frequency [GHz]", ylabel="Brightness Temperature [K]", grid=True)

Also, it is possible to execute a combination of absorption models. The following example use :code:`R19SD` model for :math:`O_2` and
:code:`R16` for :math:`H_2O`: to compute upwelling brightness temperature using emissivity surface.

.. plot::
    :include-source: true

    from pyrtlib.tb_spectrum import TbCloudRTE
    from pyrtlib.absorption_model import H2OAbsModel
    from pyrtlib.climatology import AtmosphericProfiles as atmp
    from pyrtlib.utils import ppmv2gkg, mr2rh

    z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100

    ang = np.array([90.])
    frq = np.arange(20, 201, 1)

    rte = TbCloudRTE(z, p, t, rh, frq, ang)
    rte.emissivity = 0.9
    rte.init_absmdl('R19SD')
    H2OAbsModel.model = 'R16'
    H2OAbsModel.set_ll()
    df = rte.execute()
    df = df.set_index(frq)
    df.tbtotal.plot(figsize=(12,8), xlabel="Frequency [GHz]", ylabel="Brightness Temperature [K]", grid=True)


Standard Atmospheric Profiles
=============================

Atmospheric constituent profiles (0-120km) (suplimented with other data) [ANDERSON]_
This file was partly copied from FASCOD2 routine MLATMB 10/11/87
                                                                
The file contains 6 model profiles: 

* Model 1. Tropical                                              
* Model 2. Midlatitude Summer                                    
* Model 3. Midlatitude Winter                                    
* Model 4. Subarctic Summer                                      
* Model 5. Subarctic Winter                                      
* Model 6. U.S. Standard 
  
Each of these profile contains data at 50 atmospheric levels:  
Altitude (km), Pressure (mb), Density (cm-3), Molec. densities (ppmv):
1(:math:`H_2O`),  2(:math:`CO_2`),  3(:math:`O_3`), 4(:math:`N_2O`),   5(:math:`CO`),    6(:math:`CH_4`),   7(:math:`O_2`)
Plus suplimental profiles where available.

.. autosummary::
    :toctree: generated/

    pyrtlib.climatology.AtmosphericProfiles
    pyrtlib.climatology.ProfileExtrapolation


Example:

.. code-block:: python

    from pyrtlib.climatology import AtmosphericProfiles as atmp

    z, p, d, t, md = atmp.gl_atm(atmp.TROPICAL)
    # index of available profiles
    atmp.atm_profiles()
    {0: 'Tropical',
     1: 'Midlatitude Summer',
     2: 'Midlatitude Winter',
     3: 'Subarctic Summer',
     4: 'Subarctic Winter',
     5: 'US Standard'}


Radiative Transfer Equation
===========================

RTE functions called from :py:class:`pyrtlib.rt_equation.RTEquation`:

* :code:`bright` = compute temperature for the modified Planck radiance 
* :code:`cloudy_absorption`   = computes cloud (liquid and ice) absorption profiles
* :code:`cloud_integrated_density`   = integrates cloud water density of path ds (linear) 
* :code:`cloud_radiating_temperature`   = computes mean radiating temperature of a cloud 
* :code:`clearsky_absorption`   = computes clear-sky (:math:`H_2O` and :math:`O_2`) absorption profiles
* :code:`exponential_integration`   = integrates (ln) absorption over profile layers
* :code:`planck`   = computes modified planck radiance and related quantities
* :code:`ray_tracing`  = computes refracted path length between profile levels
* :code:`refractivity`  = computes vapor pressure and refractivity profiles
* :code:`vapor`    = computes vapor pressure and vapor density 


.. autosummary::
    :toctree: generated/

    pyrtlib.rt_equation.RTEquation


Absorption Models
=================

Computes absorption coefficient in atmosphere due to water vapor (:math:`H_2O`), oxygen in air (:math:`O_2`), ozone in air (:math:`O_3`), suspended cloud liquid water droplets and 
collision-induced power absorption coefficient (neper/km) in air ("dry continuum", mostly due to :math:`N_2`-:math:`N_2`, but also contributions from :math:`O_2`-:math:`N_2` and :math:`O_2`-:math:`O_2`)

.. autosummary::
    :toctree: generated/

    pyrtlib.absorption_model.AbsModel
    pyrtlib.absorption_model.H2OAbsModel
    pyrtlib.absorption_model.O2AbsModel
    pyrtlib.absorption_model.O3AbsModel
    pyrtlib.absorption_model.N2AbsModel
    pyrtlib.absorption_model.LiqAbsModel

To get all implemented models use the following code:

.. code-block:: python

    from pyrtlib.absorption_model import AbsModel

    AbsModel.implemented_models()
    ['R98',
     'R03',
     'R16',
     'R17',
     'R19',
     'R19SD',
     'R20',
     'R20SD',
     'R21SD',
     'R22SD']

Utility Functions
=================

The utils module contains funtions of general utility used in multiple places throughout *pyrtlib*.
 
.. autosummary::
    :toctree: generated/
    :template: custom-module-template.rst

    pyrtlib.utils


Uncertainty
===========

This module has some tool to compute the absorption model sensitivity to the uncertainty of spectroscopic parameters, 
with the purpose of identifying the most significant contributions to the total uncertainty of modeled upwelling/downwelling
brightness temperture.

.. autosummary::
    :toctree: generated/

    pyrtlib.uncertainty.AbsModUncertainty
    pyrtlib.uncertainty.SpectroscopicParameter

API Web Services
================
Observations dataset web services which may be used in pyrtlib. 
Available datasets are the Wyoming Upper Air Archive (University of Wyoming), NCEI’s Integrated Radiosonde Archive version 2 (IGRA2) or the 
ERA5 Reanalysis model data (Copernicus Climate Change Service). See examples to get started to use these services.

.. note::
    Parts of the code have been reused from the `Siphon <https://github.com/Unidata/siphon>`_ library.

.. autosummary::
    :toctree: generated/

    pyrtlib.apiwebservices.WyomingUpperAir
    pyrtlib.apiwebservices.IGRAUpperAir
    pyrtlib.apiwebservices.ERA5Reanalysis
    
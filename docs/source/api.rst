==============
API references
==============

Main class
====================

The main class which computes brightness temperatures (Tb), mean
radiating temperature (Tmr), and integrated absorption (Tau) for 
clear or cloudy conditions,  Also returns all integrated quantities
that the original TBMODEL, Cyber Version, returned.

.. autosummary::
    :toctree: generated/

    pyrtlib.main.BTCloudRTE

Example
.......

.. code-block:: python

    from pyrtlib.main import BTCloudRTE

    rte = BTCloudRTE(z, p, t, rh, frq, ang)
    rte.init_absmdl('rose19sd')
    rte.satellite = True
    rte.emissivity = 0.6
    df = rte.execute()

Also, it is possible to execute a combination of absorption models. The following example use :code:`rose19sd` model for :math:`O_2` and
:code:`rose16` for :math:`H_2O`:

.. code-block:: python

    from pyrtlib.main import BTCloudRTE
    from pyrtlib.utils import import_lineshape
    from pyrtlib.absmodel import H2OAbsModel, O2AbsModel

    rte = BTCloudRTE(z, p, t, rh, frq, ang)
    rte.init_absmdl('rose19sd')
    H2OAbsModel.model = 'rose16'
    H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
    df = rte.execute()


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

    pyrtlib.atmp.AtmosphericProfiles


Example
.......

.. code-block:: python

    from pyrtlib.atmp import AtmosphericProfiles as atmp

    z, p, d, tk, md = atmp.gl_atm(atmp.TROPICAL)
    # index of available profiles
    atmp.atm_profiles()


Radiative Transfer Equation
===========================

RTE functions called from :py:class:`pyrtlib.rte.RTEquation`:

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

    pyrtlib.rte.RTEquation


Absorption Models
=================

Computes absorption coefficient in atmosphere due to water vapor (:math:`H_2O`), oxygen in air (:math:`O_2`), ozone in air (:math:`O_3`), suspended cloud liquid water droplets and 
collision-induced power absorption coefficient (neper/km) in air ("dry continuum", mostly due to :math:`N_2`-:math:`N_2`, but also contributions from :math:`O_2`-:math:`N_2` and :math:`O_2`-:math:`O_2`)

.. autosummary::
    :toctree: generated/
    :template: custom-class-template.rst

    pyrtlib.absmodel.AbsModel
    pyrtlib.absmodel.H2OAbsModel
    pyrtlib.absmodel.O2AbsModel
    pyrtlib.absmodel.O3AbsModel
    pyrtlib.absmodel.N2AbsModel
    pyrtlib.absmodel.LiqAbsModel


Utility Functions
=================

The utils module contains funtions of general utility used in multiple places throughout *pyrtlib*.
 
.. autosummary::
    :toctree: generated/
    :template: custom-module-template.rst

    pyrtlib.utils


Uncertainty
===========

.. autosummary::
    :toctree: generated/
    :template: custom-module-template.rst

    pyrtlib.lineshape
    pyrtlib.absmod_uncertainty


API Web Services
================
Observations dataset web services which may be used in pyrtlib. 
Available datasets are the Wyoming Upper Air Archive (University of Wyoming), NCEI’s Integrated Radiosonde Archive version 2 or the 
ERA5 Reanalysis model data (Copernicus Climate Change Service). See examples to get started to use these services.

.. autosummary::
    :toctree: generated/

    pyrtlib.apiwebservices.WyomingUpperAir
    pyrtlib.apiwebservices.IGRAUpperAir
    pyrtlib.apiwebservices.ERA5Reanalysis
    
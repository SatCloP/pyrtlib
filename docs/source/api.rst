##############
API references
##############

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

Atmospheric Profiles
====================

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
1(H2O),  2(CO2),  3(O3), 4(N2O),   5(CO),    6(CH4),   7(O2)
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

.. autosummary::
    :toctree: generated/

    pyrtlib.rte.RTEquation


Absorption Models
=================

.. autosummary::
    :toctree: generated/
    :template: custom-class-template.rst

    pyrtlib.absmodel.AbsModel
    pyrtlib.absmodel.H2OAbsModel
    pyrtlib.absmodel.O2AbsModel
    pyrtlib.absmodel.N2AbsModel
    pyrtlib.absmodel.LiqAbsModel


Utility Function
================
 
.. autosummary::
    :toctree: generated/
    :template: custom-module-template.rst   

    pyrtlib.utils


Line Shape
===========

.. autosummary::
    :toctree: generated/

    pyrtlib.lineshape

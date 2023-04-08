.. .. image:: ../../resources/logo/logo_large_new.png
..   :width: 600

.. pyrtlib documentation master file, created by sphinx-quickstart on Fri Mar 19 09:49:16 2021. 
   You can adapt this file completely to your liking, but it should at least contain the root `toctree` directive.


PyRTLib documentation
=====================

PyRTLib allows to simulate and calculate radiometric parameters and estimting propogation parameters using as input meteorological data.
Some meteorological dataset are built-in in PyRTLib which can be download and used directly in PyRTLib. It considers atmospheric profiles from both radiosounding observations (RAOB) and model reanalysis (ERA5).
RAOB profiles come from Wyoming Upper Air Archive (University of Wyoming) and NCEI’s Integrated Radiosonde Archive version 2 by the National Climatic Data Center (NCDC) of the National Oceanic and Atmospheric Administration (NOAA).

PyRTLib also allows to quantify absorption model uncertainty due to uncertainty in the underlying spectroscopic parameters. [Cimini-2018]_
The approach is applied to a widely used microwave absorption model [Rosenkranz-2017]_, on which PyRTLib is based, and radiative transfer calculations at any frequencies range, 
which are commonly exploited for atmospheric sounding by microwave radiometer (MWR).  

.. panels::
   :card: + intro-card text-center
   :column: col-lg-6 col-md-6 col-sm-6 col-xs-12 d-flex
   
   ---
   :img-top: _static/gallery_panels.svg
   
   Gallery example

   +++

   .. link-button:: examples/index
      :type: ref
      :text: Go To Reference
      :classes: btn-outline-primary btn-block
   
   ---
   :img-top: _static/gallery_panels.svg

   Jupyter notebook example

   +++

   .. link-button:: notebook/index
      :type: ref
      :text: Go To Reference
      :classes: btn-block btn-outline-primary 

.. pyrtlib is a python tool that provides a set of calsses and methods for simulating ........

The source code for pyrtlib python package is hosted on `github
<https://github.com/slarosa/pyrtlib>`_.

.. cssclass:: image-pyrtlib
   
   .. image:: ../../resources/spectrum_r22.jpeg
      :width: 600
   .. image:: ../../resources/r98_r22.jpeg
      :width: 600

Example:
--------
Atmospheric profile definition:

>>> z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

Units conversion:

>>> gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)

Relative humidity of H2O (water vapor)

>>> rh = mr2rh(p, t, gkg)[0] / 100

Deifinition of angles and frequencies:

>>> ang = np.array([90.])
>>> frq = np.arange(20, 201, 1)

Initialize parameters for main execution:

>>> rte = TbCloudRTE(z, p, t, rh, frq, ang)

Set absorption model:

>>> rte.init_absmdl('rose16')

Execute model by computing upwelling radiances:

>>> df = rte.execute()
>>> df.tbtotal
0      297.391838
1      296.186240
2      294.748245
3      294.953483
4      296.027799
         ...
176    275.997899
177    276.611319
178    277.129218
179    277.566840
180    277.936645
Name: tbtotal, Length: 181, dtype: float64

Preview of the output dataframe (see :py:meth:`pyrtlib.tb_spectrum.TbCloudRTE.execute` for more info):

.. list-table::
   :widths: 25 25 25 25 25 25 25 25
   :header-rows: 1

   * - tbtotal
     - tbatm
     - tmr
     - tmrcld
     - tauwet
     - taudry
     - tauliq
     - tauice
   * - 297.391838 
     - 0.0 
     - 281.191287
     - 0.0
     - 0.120341 
     - 0.012855
     - 0.0
     - 0.0
   * - 296.186240    
     - 0.0      
     - 280.517137    
     - 0.0       
     - 0.188802    
     - 0.013524    
     - 0.0       
     - 0.0
   * - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
     - ...
   * - 277.936645    
     - 0.0      
     - 276.874793    
     - 0.0       
     - 3.041382    
     - 0.026476    
     - 0.0      
     -  0.0 
  
.. toctree::
   :maxdepth: 2
   :caption: Installation:

   installation

.. toctree::
   :maxdepth: 2
   :caption: API References:
   
   api

.. toctree::
   :maxdepth: 2
   :caption: Notebook:
   
   notebook/index

.. toctree::
   :maxdepth: 2
   :caption: Examples:
   
   examples/index

.. toctree::
   :maxdepth: 1
   :caption: References:
   
   references

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |build-docs-action| image:: https://github.com/slarosa/pyrtlib/workflows/build-docs-action/badge.svg
   :target: https://github.com/slarosa/pyrtlib/actions/workflows/build_docs.yml

.. |run-python-tests| image:: https://github.com/slarosa/pyrtlib/workflows/run-python-tests/badge.svg
   :target: https://github.com/slarosa/pyrtlib/actions/workflows/ci.yml

.. |license| image:: https://img.shields.io/github/license/slarosa/pyrtlib.svg
   :target: https://github.com/slarosa/pyrtlib/blob/main/LICENSE.md

.. |GitHub commit| image:: https://img.shields.io/github/last-commit/slarosa/pyrtlib
   :target: https://github.com/slarosa/pyrtlib/commits/main

.. |codecov| image:: https://codecov.io/gh/slarosa/pyrtlib/branch/main/graph/badge.svg?token=7DV4B4U1OZ
   :target: https://codecov.io/gh/slarosa/pyrtlib
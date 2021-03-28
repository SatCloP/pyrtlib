.. image:: ../../resources/logo/logo_large.png
  :width: 600

.. pyrtlib documentation master file, created by
   sphinx-quickstart on Fri Mar 19 09:49:16 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyrtlib's documentation!
====================================================

pyrtlib is a python tool that provides a set of
functions and classes for simulating ........

The source code for pyrtlib python package is hosted on `github
<https://github.com/slarosa/pyrtlib>`_.

.. image:: ../../resources/spectrum.png

.. code-block:: python
   :linenos:

   import numpy as np
   from pyrtlib.rte import RTEquation
   from pyrtlib.absmodel import H2OAbsModel, O2AbsModel
   from pyrtlib.linelist import h2o_linelist, o2_linelist
   from pyrtlib.atmp import AtmosphericProfiles as atmp
   from pyrtlib.utils import ppmv2gkg, mr2rh

   z, p, d, tk, md = atmp.gl_atm(atmp.TROPICAL)
   frq = np.arange(20, 201, 1)
   ice = 0
   gkg = ppmv2gkg(md[:, atmp.H20], atmp.H20)
   rh = mr2rh(p, tk, gkg)[0] / 100

   e, rho = RTEquation.vapor(tk, rh, ice)

   H2OAbsModel.model = 'rose19sd'
   H2OAbsModel.h2oll = h2o_linelist()
   O2AbsModel.model = 'rose19sd'
   O2AbsModel.o2ll = o2_linelist()
   for i in range(0, len(z)):
      v = 300.0 / tk[i]
      ekpa = e[i] / 10.0
      pdrykpa = p[i] / 10.0 - ekpa
      for j in range(0, len(frq)):
         _, _ = H2OAbsModel().h2o_rosen19_sd(pdrykpa, v, ekpa, frq[j])
         # _, _ = O2AbsModel().o2abs_rosen18(pdrykpa, v, ekpa, frq[j])

   # for j in range(0, len(frq)):
   #     awet, adry = RTEquation.clearsky_absorption(p, tk, e, frq[j])
   #     npp, ncpp = H2OAbsModel().h2o_rosen19_sd(pdrykpa, v, ekpa, frq)

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

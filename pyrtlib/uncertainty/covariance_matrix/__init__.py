# This module contains uncertainty covariance matrix for oxygen and water vapor parameters
# as estimated in Cimini et al. ACP 2018 (https://doi.org/10.5194/acp-18-15231-2018) and
# updated in Cimini et al. GMD 2019 (https://doi.org/10.5194/gmd-12-1833-2019).
# The covariance matrices were derived with the routines by Phil Rosenkranz (https://doi.org/10.21982/M81013),
# as downloaded on 15 June 2017, available at http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html
#
# This version is VALID for ground-based RT calculations in the range 20-150 GHz.
#
# The list of 112 parameters with units is as follows:
#   1      : O2 S(300) [%]^2
#   2      : O2 n_a [adim]^2
#   3      : O2 gamma_0(300) [GHz/bar]^2
#   4- 37  : O2 gamma_a(300) [GHz/bar]^2 for N=1-,1+,3-,...,33-,33+ lines
#  38- 71  : O2 y(300) [1/bar]^2 for N=1-,1+,3-,...,33-,33+ lines
#  72-105  : O2 v [1/bar]^2 for N=1-,1+,3-,...,33-,33+ lines
# 106      : H2O C_f(300) [km-1 mb-2 GHz-2]^2
# 107      : H2O C_s(300) [km-1 mb-2 GHz-2]^2
# 108      : H2O gamma_a(296) at 22.2 GHz [GHz/bar]^2
# 109      : H2O S(296) at 22.2 GHz [Hz*cm2]^2
# 110      : H2O n_Cf [adim]^2
# 111      : H2O R at 22.2 GHz [adim]^2
# 112      : H2O n_Cs [adim]^2

import os

import numpy as np
from netCDF4 import Dataset

PATH = os.path.dirname(os.path.abspath(__file__))

nc = Dataset(os.path.join(PATH, "R17", "Cov_parameters_Cimini_et_al_2018_V1.1.nc"))
R17_111 = np.asarray(nc.variables['Cov_p'][:])

nc = Dataset(os.path.join(PATH, "R17", "Cov_parameters_Cimini_et_al_2019_V2.0.nc"))
R17_112 = np.asarray(nc.variables['Cov_p'][:])

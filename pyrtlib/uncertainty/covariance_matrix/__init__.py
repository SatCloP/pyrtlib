# This file contains uncertainty covariance matrix for oxygen and water vapor parameters
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

nc = Dataset(os.path.join(
    PATH, "R17", "Cov_parameters_Cimini_et_al_2018_V1.1.nc"))
R17_111 = np.asarray(nc.variables['Cov_p'][:])
nc.close()

nc = Dataset(os.path.join(
    PATH, "R17", "Cov_parameters_Cimini_et_al_2019_V2.0.nc"))
R17_112 = np.asarray(nc.variables['Cov_p'][:])
nc.close()

# This file contains uncertainty covariance matrix for oxygen and water vapor parameters
# as estimated in Gallucci et al. ACPD 2023 (https://doi.org/10.5194/egusphere-2023-3160).
# The reference code is Rosenkranz 2019 - version of 2019-02-01,
# available at: http://cetemps.aquila.infn.it/mwrnet/lblmrt_ns.html,
# Note that this is an extension to the one developed in Cimini et al. ACP 2018,
# https://doi.org/10.5194/acp-18-15231-2018-supplement,
# which is also valid for upwelling Brightness Temperatures,
# Radiative Transfer calculations in the range 16-700 GHz.
#
# The list of 135 parameters with units is as follows:
# 1        : O2  S(300) [%]
# 2        : O2  n_a [adim]
# 3        : O2  gamma_0(300) [GHz/bar]
# 4- 37    : O2  gamma_a(300) [GHz/bar] for N=1-,1+,3-,...,33-,33+ lines
# 38- 71   : O2  y(300) [1/bar] for N=1-,1+,3-,...,33-,33+ lines
# 72-105   : O2  v [1/bar] for N=1-,1+,3-,...,33-,33+ lines
# 106      : H2O C_f(300) [km-1 mb-2 GHz-2]
# 107      : H2O C_s(300) [km-1 mb-2 GHz-2]
# 108      : H2O n_Cf [adim]
# 109      : H2O n_Cs [adim]
# 110      : H2O S 022GHz (line 01)    [Hz*cm2]
# 111      : H2O S 183GHz (line 02)    [Hz*cm2]
# 112      : H2O S 325GHz (line 04)    [Hz*cm2]
# 113      : H2O S 448GHz (line 08)    [Hz*cm2]
# 114      : H2O S 556GHz (line 12)    [Hz*cm2]
# 115      : H2O S 752GHz (line 15)    [Hz*cm2]
# 116      : H2O gamma_a 022GHz line01 [GHz/bar]
# 117      : H2O gamma_a 183GHz line02 [GHz/bar]
# 118      : H2O gamma_a 325GHz line04 [GHz/bar]
# 119      : H2O gamma_a 448GHz line08 [GHz/bar]
# 120      : H2O gamma_a 556GHz line12 [GHz/bar]
# 121      : H2O gamma_a 752GHz line15 [GHz/bar]
# 122      : H2O n_a 022GHz (line 01)  [adim]
# 123      : H2O n_a 183GHz (line 02)  [adim]
# 124      : H2O n_a 325GHz (line 04)  [adim]
# 125      : H2O n_a 448GHz (line 08)  [adim]
# 126      : H2O n_a 556GHz (line 12)  [adim]
# 127      : H2O n_a 752GHz (line 15)  [adim]
# 128      : H2O S 380GHz (line 05)    [Hz*cm2]
# 129      : H2O S 474GHz (line 10)    [Hz*cm2]
# 130      : H2O S 620GHz (line 13)    [Hz*cm2]
# 131      : H2O gamma_a 620GHz line13 [GHz/bar]
# 132      : O2  gamma_a(300) 234 GHz   [GHz/bar]
# 133      : O2  gamma_a(300) 368 GHz   [GHz/bar]
# 134      : O2  gamma_a(300) 424 GHz   [GHz/bar]
# 135      : O2  gamma_a(300) 487 GHz   [GHz/bar]

PATH = os.path.dirname(os.path.abspath(__file__))

nc = Dataset(os.path.join(
    PATH, "R19", "Cov_parameters_Gallucci_et_al_ACPD_2023_V1.0.nc"))
R19 = np.asarray(nc.variables['Cov_p'][:])
nc.close()

# -*- coding: utf-8 -*-

"""REFERENCES FOR MEASUREMENTS (freq in GHz, S in Hz*cm^2, W's,D's in GHZ/bar, X's & A dimensionless) updated Dec. 30, 2020.

    References
    ----------
    .. [1] M. Koshelev et al., JQSRT v.205, pp. 51-58 (2018)
    .. [2] V. Payne et al.,IEEE Trans. Geosci. Rem. Sens. v.46, pp.3601-3617 (2008)
    .. [3] G. Golubiatnikov, J. MOLEC. SPEC. vol. 230, pp.196-198 (2005)
    .. [4] M. Koshelev et al., J. Molec. Spec. v.241, pp.101-108 (2007)
    .. [5] J.-M. Colmont et al.,J. Molec. Spec. v.193, pp.233-243 (1999)
    .. [6] M. Tretyakov et al, JQSRT v.114 pp.109-121 (2013)
    .. [7] G. Golubiatnikov et al., JQSRT v.109, pp.1828-1833 (2008)
    .. [8] V. Podobedov et al., JQSRT v.87, pp. 377-385 (2004)
    .. [9] M. Koshelev, JQSRT v.112, pp.550-552 (2011)
    .. [10] M. Tretyakov, JQSRT v.328, pp.7-26 (2016)
    .. [11] D. Turner et al., IEEE Trans. Geosci. Rem. Sens. v.47 pp.3326-37 (2009), re-adjusted for new line par. Aug.22, 2022.
    .. [11] M. Koshelev et al. JQSRT doi:10.1016/j.jqsrt.2020.107472

    Other parameters from HITRAN2020.
"""
# code to import coefficient from .asc file
# import pandas as pd
# import numpy as np

# a = pd.read_table("h2o_sdlist.asc", sep=',',header=None, skiprows=1, nrows=19, usecols=range(0, 20))
# a = np.loadtxt("h2o_sdlist.asc", delimiter=',', skiprows=1, usecols=range(0, 20), max_rows=20)
# np.savetxt("/Users/slarosa/dev/pyrtlib/pyrtlib/lineshape/h2o_list_rose22.txt", a)

import numpy as np
import os

PATH = os.path.dirname(os.path.abspath(__file__))
mtx = np.loadtxt(os.path.join(PATH, "h2o_list_r22.txt"))

# import numpy as np

# # molecule freq,GHz  S(296K)    B     W0air XWair W0self XWself  Dair  XDair  Dself XDself  Aair
# # Aself W2air W2self Refs. FL(i)     S1(i)    B2(i)   Wair  X(i)  Wself  Xs(i)   Sair  Xh(I)  Sself
# # Xhs(I) Aair Aself  W2   W2S (variable names in the code) Refs.
# blk = np.nan

# mtx = np.vstack(
#     [[11, 22.23508, 1.335e-14, 2.172, 2.74, 0.76, 13.63, 1.2, -0.033, 2.6, 0.814, blk, 0.0, 0.0, 0.435, blk, 1.91, blk,
#       0.0, 0.0],
#      [11, 183.310087, 2.319e-12, 0.677, 3.033, 0.62, 15.01, 0.82, -0.074, 1.8, 0.136, 0.98, 0.0, 12.6, 0.407, 0.412,
#       1.46, 0.571, -0.016, 0.16],
#      [11, 321.22563, 7.657e-14, 6.262, 2.426, 0.73, 10.65, 0.54, -0.143, blk, 0.278, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 325.152888, 2.721e-12, 1.561, 2.847, 0.64, 13.95, 0.74, -0.013, blk, 1.325, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 380.197353, 2.477e-11, 1.062, 2.868, 0.54, 14.4, 0.89, -0.074, blk, 0.24, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 439.150807, 2.137e-12, 3.643, 2.055, 0.69, 9.06, 0.52, 0.051, blk, 0.165, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 443.018343, 4.44e-13, 5.116, 1.819, 0.7, 7.96, 0.5, 0.14, blk, -0.229, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 448.001085, 2.588e-11, 1.424, 2.612, 0.7, 13.01, 0.67, -0.116, blk, -0.615, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 470.888999, 8.196e-13, 3.645, 2.169, 0.73, 9.7, 0.65, 0.061, blk, -0.465, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 474.689092, 3.268e-12, 2.411, 2.366, 0.71, 11.24, 0.64, -0.027, blk, -0.72, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 488.490108, 6.628e-13, 2.89, 2.616, 0.75, 13.58, 0.72, -0.065, blk, -0.36, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 556.935985, 1.57e-09, 0.161, 3.115, 0.75, 14.24, 1.0, 0.187, blk, -1.693, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 620.700807, 1.7e-11, 2.423, 2.468, 0.79, 11.94, 0.75, 0.0, blk, 0.687, 0.92, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 658.006072, 9.033e-13, 7.921, 3.154, 0.73, 13.84, 1.0, 0.176, blk, -1.496, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk],
#      [11, 752.033113, 1.035e-09, 0.402, 3.114, 0.77, 13.58, 0.84, 0.162, blk, -0.878, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 859.965608, 0.5705e-12, 8.163, 3.121, 0.76, 14.08, 0.76, 0.005, blk, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 916.171582, 4.275e-11, 1.461, 2.695, 0.79, 13.55, 0.48, -.001, blk, 0.521, 0.47, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 970.315045, 0.4806E-10, 1.944, 2.574, 0.70, 25.95, 0.70, -.003, blk, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 987.926803, 0.7528E-09, 0.261, 2.976, 0.75, 14.35, 0.75, -.002, blk, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk],
#      [11, 1097.36487, 0.4890E-08, 0.754, 3.095, 0.75, 15.27, 0.75, 0.002, blk, blk, blk, blk, blk, blk, blk, blk,
#       blk, blk, blk]
#      ])
# # continuum terms
ctr = np.array([300.0, 5.919e-10, 3.0, 1.416e-08, 7.5])

# # below is from abh2o_sd.f
# # read line parameters; units: ghz, hz*cm^2, mhz/mb
reftline = 296.0
fl = mtx[:, 1]
s1 = mtx[:, 2]
b2 = mtx[:, 3]
w0 = mtx[:, 4] / 1000.0
x = mtx[:, 5]
w0s = mtx[:, 6] / 1000.0
xs = mtx[:, 7]
sh = mtx[:, 8] / 1000.0
xh = mtx[:, 9]
shs = mtx[:, 10] / 1000.0
xhs = mtx[:, 11]
aair = mtx[:, 12]
aself = mtx[:, 13]
w2 = mtx[:, 14] / 1000.0
xw2 = mtx[:, 15]
w2s = mtx[:, 16] / 1000.0
xw2s = mtx[:, 17]
d2 = mtx[:, 18] / 1000.0
d2s = mtx[:, 19] / 1000.0

# # replace non-existing shifting parameters with broadening parameters
# indx = np.where(np.isnan(xh))
# xh[indx] = x[indx]
# indx = np.where(np.isnan(xhs))
# xhs[indx] = xs[indx]
# indx = np.where(np.isnan(xw2))
# xw2[indx] = x[indx]
# indx = np.where(np.isnan(xw2s))
# xw2s[indx] = xs[indx]
# # replace non-existing aair aself parameters with zero (to be agreed with phil)
# indx = np.where(np.isnan(aair))
# aair[indx] = 0
# indx = np.where(np.isnan(aself))
# aself[indx] = 0
# indx = np.where(np.isnan(w2))
# w2[indx] = 0
# indx = np.where(np.isnan(w2s))
# w2s[indx] = 0
# indx = np.where(np.isnan(d2))
# d2[indx] = 0
# indx = np.where(np.isnan(d2s))
# d2s[indx] = 0
# indx = np.where(np.isnan(sh))
# sh[indx] = 0
# indx = np.where(np.isnan(shs))
# shs[indx] = 0
# # read continuum parameters; units: kelvin, 1/(km*mb^2*ghz^2)
reftcon = ctr[0]
cf = ctr[1]
xcf = ctr[2]
cs = ctr[3]
xcs = ctr[4]


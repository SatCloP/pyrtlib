import numpy as np

# molecule freq,GHz  S(296K)    B     W0air XWair W0self XWself  Dair  XDair  Dself XDself  Aair
# Aself W2air W2self Refs. FL(i)     S1(i)    B2(i)   Wair  X(i)  Wself  Xs(i)   Sair  Xh(I)  Sself
# Xhs(I) Aair Aself  W2   W2S (variable names in the code) Refs.
blk = np.nan

mtx = np.vstack([[11, 22.23508, 1.335e-14, 2.172, 2.699, 0.76, 13.29, 1.2, - 0.033, 2.6, 0.814, blk, blk, blk],
                 [11, 183.310087, 2.319e-12, 0.677, 2.959, 0.63, 14.81, 0.83, - 0.072, 1.8, 0.108, 1.25, 0.0, 17.3],
                 [11, 321.22563, 7.657e-14, 6.262, 2.426, 0.73, 10.65, 0.54, - 0.143, blk, 0.278, blk, blk, blk],
                 [11, 325.152888, 2.721e-12, 1.561, 2.847, 0.64, 13.95, 0.74, - 0.013, blk, 1.325, blk, blk, blk],
                 [11, 380.197353, 2.477e-11, 1.062, 2.868, 0.54, 14.4, 0.89, - 0.074, blk, 0.24, blk, blk, blk],
                 [11, 439.150807, 2.137e-12, 3.643, 2.055, 0.69, 9.06, 0.52, 0.051, blk, 0.165, blk, blk, blk],
                 [11, 443.018343, 4.44e-13, 5.116, 1.819, 0.7, 7.96, 0.5, 0.14, blk, - 0.229, blk, blk, blk],
                 [11, 448.001085, 2.588e-11, 1.424, 2.612, 0.7, 13.01, 0.67, - 0.116, blk, - 0.615, blk, blk, blk],
                 [11, 470.888999, 8.196e-13, 3.645, 2.169, 0.73, 9.7, 0.65, 0.061, blk, - 0.465, blk, blk, blk],
                 [11, 474.689092, 3.268e-12, 2.411, 2.366, 0.71, 11.24, 0.64, - 0.027, blk, - 0.72, blk, blk, blk],
                 [11, 488.490108, 6.628e-13, 2.89, 2.616, 0.75, 13.58, 0.72, - 0.065, blk, - 0.36, blk, blk, blk],
                 [11, 556.935985, 1.57e-09, 0.161, 3.115, 0.75, 14.24, 1.0, 0.187, blk, - 1.693, blk, blk, blk],
                 [11, 620.700807, 1.7e-11, 2.423, 2.468, 0.79, 11.94, 0.75, 0.0, blk, 0.687, 0.92, blk, blk],
                 [11, 658.006072, 9.033e-13, 7.921, 3.154, 0.73, 13.84, 1.0, 0.176, blk, - 1.496, blk, blk, blk],
                 [11, 752.033113, 1.035e-09, 0.402, 3.114, 0.77, 13.58, 0.84, 0.162, blk, - 0.878, blk, blk, blk],
                 [11, 916.171582, 4.275e-11, 1.461, 2.695, 0.79, 13.55, 0.48, 0.0, blk, 0.521, 0.47, blk, blk]])
# continuum terms
ctr = np.array([300.0, 5.954e-10, 3.0, 1.42e-08, 7.5])

# below is from abh2o_sd.f
# read line parameters; units: ghz, hz*cm^2, mhz/mb
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

# replace non-existing shifting parameters with broadening parameters
indx = np.where(np.isnan(xh))
xh[indx] = x[indx]
indx = np.where(np.isnan(xhs))
xhs[indx] = xs[indx]
# replace non-existing aair aself parameters with zero (to be agreed with phil)
indx = np.where(np.isnan(aair))
aair[indx] = 0
indx = np.where(np.isnan(aself))
aself[indx] = 0
# read continuum parameters; units: kelvin, 1/(km*mb^2*ghz^2)
reftcon = ctr[0]
cf = ctr[1]
xcf = ctr[2]
cs = ctr[3]
xcs = ctr[4]

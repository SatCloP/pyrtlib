import numpy as np
import os

PATH = os.path.dirname(os.path.abspath(__file__))
mtx = np.loadtxt(os.path.join(PATH, "h2o_list_r21.txt"))
# continuum terms
ctr = np.array([300.0, 5.919e-10, 3.0, 1.416e-08, 7.5])

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
w2 = mtx[:, 14] / 1000.0
xw2 = mtx[:, 15]
w2s = mtx[:, 16] / 1000.0
xw2s = mtx[:, 17]
d2 = mtx[:, 18] / 1000.0
d2s = mtx[:, 19] / 1000.0

reftcon = ctr[0]
cf = ctr[1]
xcf = ctr[2]
cs = ctr[3]
xcs = ctr[4]

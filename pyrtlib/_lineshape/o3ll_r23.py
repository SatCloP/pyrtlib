import numpy as np
import os

# a = pd.read_table("o3_list.asc",header=None, skiprows=1, nrows=464, usecols=range(0, 7), delim_whitespace=True)
# a = np.loadtxt("h2o_sdlist.asc", delimiter=',', skiprows=1, usecols=range(0, 20), max_rows=20)
# np.savetxt("o3_list_r23.txt", a)

PATH = os.path.dirname(os.path.abspath(__file__))
os.path.join(PATH, "o3_list_r23.txt")
mtx = np.loadtxt(os.path.join(PATH, "o3_list_r23.txt"))

reftline = 296
# line intensities from hitran include isotopomer abundance
# read line parameters; units: ghz, hz*cm^2, mhz/mb
fl = mtx[:, 1]
s1 = mtx[:, 2]
b = mtx[:, 3]
w = mtx[:, 4] / 1000.0
x = mtx[:, 5]
sr = mtx[:, 6]

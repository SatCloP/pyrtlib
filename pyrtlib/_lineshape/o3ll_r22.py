import numpy as np
import os

PATH = os.path.dirname(os.path.abspath(__file__))
os.path.join(PATH, "o3_list_r22.txt")
mtx = np.loadtxt(os.path.join(PATH, "o3_list_r22.txt"))

reftline = 296
# line intensities from hitran include isotopomer abundance
# read line parameters; units: ghz, hz*cm^2, mhz/mb
fl = mtx[:, 1]
s1 = mtx[:, 2]
b = mtx[:, 3]
w = mtx[:, 4] / 1000.0
x = mtx[:, 5]
sr = mtx[:, 6]

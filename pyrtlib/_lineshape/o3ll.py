# -*- coding: utf-8 -*-

# nc = Dataset("/Users/slarosa/Downloads/o3_lineshape.nc", mode='w')

# for m in O3AbsModel.implemented_models():
#     if m in ['R22SD', 'R18']:
#         if m == 'R22SD':
#             m = 'R22'
#         grp = nc.createGroup(m)
#         O3AbsModel.model = m
#         O3AbsModel.set_ll()
#         grp.createDimension('rows', O3AbsModel.o3ll.mtx.shape[0])
#         grp.createDimension('cols', O3AbsModel.o3ll.mtx.shape[1])
#         mtx = grp.createVariable('mtx', 'f8', ('rows', 'cols', ))
#         mtx[:] = O3AbsModel.o3ll.mtx
#         reftline = grp.createVariable('reftline', 'f4', ())
#         reftline[:] = O3AbsModel.o3ll.reftline

# nc.close()

import os
from netCDF4 import Dataset

from pyrtlib.absorption_model import O3AbsModel

PATH = os.path.dirname(os.path.abspath(__file__))
nc = Dataset(os.path.join(PATH, "o3_lineshape.nc"), mode='r')

d = nc.groups[O3AbsModel.model]
mtx = d.variables['mtx'][:].data
reftline = d.variables['reftline'][:].data.item()

fl = mtx[:, 1]
s1 = mtx[:, 2]
b = mtx[:, 3]
w = mtx[:, 4] / 1000.0
x = mtx[:, 5]
sr = mtx[:, 6]

nc.close()

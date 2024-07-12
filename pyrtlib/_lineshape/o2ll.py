# -*- coding: utf-8 -*-

# nc = Dataset("/Users/slarosa/Downloads/o2_lineshape.nc", mode='w')

# for m in O2AbsModel.implemented_models():
#     if m in ['R21SD']:
#         continue
#     O2AbsModel.model = m
#     try:
#         O2AbsModel.set_ll()
#     except:
#         m = m.replace("SD", "")
#         O2AbsModel.model = m.replace("SD", "")
#         O2AbsModel.set_ll()
#     grp = nc.createGroup(m)
#     vl = [item for item in dir(O2AbsModel.o2ll) if not item.startswith("__")]
#     for v in vl:
#         if v != 'np':
#             if v in ['wb300', 'x']:
#                 mtx = grp.createVariable(v, 'f4', ())
#             else:
#                 grp.createDimension(f'rows_{v}', getattr(O2AbsModel.o2ll, v).shape[0])
#                 mtx = grp.createVariable(v, 'f8', (f'rows_{v}',))
#             mtx[:] = getattr(O2AbsModel.o2ll, v)

# nc.close()

import os
from netCDF4 import Dataset

from pyrtlib.absorption_model import O2AbsModel

PATH = os.path.dirname(os.path.abspath(__file__))
nc = Dataset(os.path.join(PATH, "o2_lineshape.nc"), mode='r')

d = nc.groups[O2AbsModel.model]

f = d.variables['f'][:].data
s300 = d.variables['s300'][:].data
be = d.variables['be'][:].data
wb300 = d.variables['wb300'][:].data.item()
x = d.variables['x'][:].data.item()
w300 = d.variables['w300'][:].data
if O2AbsModel.model in ['R98', 'R03', 'R17', 'R18', 'R19', 'R19SD']:
    v = d.variables['v'][:].data
if O2AbsModel.model in ['R98', 'R03', 'R17', 'R18', 'R19', 'R19SD', 'R23', 'R24']:
    y300 = d.variables['y300'][:].data
if O2AbsModel.model not in ['R98', 'R03', 'R17', 'R18', 'R19', 'R19SD', 'R23', 'R24']:
    y0 = d.variables['y0'][:].data
if O2AbsModel.model in ['R16', 'R20', 'R20SD', 'R21', 'R22', 'R23', 'R24']:
    y1 = d.variables['y1'][:].data
    g0 = d.variables['g0'][:].data
    g1 = d.variables['g1'][:].data
    dnu0 = d.variables['dnu0'][:].data
    dnu1 = d.variables['dnu1'][:].data

nc.close()

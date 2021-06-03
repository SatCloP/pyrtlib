"""
Performing Brightness Temperature calculation form satellite with Ozone
=======================================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.main.BTCloudRTE` method to calculate brightness temperature from satellite with ozone.

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np

from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.absmodel import H2OAbsModel, O3AbsModel
from pyrtlib.utils import ppmv2gkg, mr2rh, ppmv_to_moleculesm3, import_lineshape

atm = ['Tropical',
       'Midlatitude Summer',
       'Midlatitude Winter',
       'Subarctic Summer',
       'Subarctic Winter',
       'U.S. Standard']

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('${T_B}$ [K]')

z, p, d, t, md = atmp.gl_atm(5) # 'U.S. Standard'

o3n_ppmv = md[:, atmp.O3]
o3n = np.zeros(z.shape)
for k in range(0, len(z)):
    o3n[k] = ppmv_to_moleculesm3(o3n_ppmv[k], p[k] * 100.0, t[k])
# o3n = np.interp(z,o3n,z)
# o3n_matlab = np.loadtxt("/Users/slarosa/Downloads/o3.csv", delimiter=',')
# assert_allclose(o3n, o3n_matlab, atol=0)

gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
rh = mr2rh(p, t, gkg)[0] / 100

ang = np.array([90.])
frq = np.arange(20, 201, 1)

rte = BTCloudRTE(z, p, t, rh, frq, ang, o3n)
rte.init_absmdl('rose20')
H2OAbsModel.model = 'rose21sd'
H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
O3AbsModel.model = 'rose18'
O3AbsModel.o3ll = import_lineshape('o3ll_{}'.format(O3AbsModel.model))
df = rte.execute()

df = df.set_index(frq)
df.tbtotal.plot(ax=ax, linewidth=1, label='{} - {}'.format(atm[5], 'rose21sd'))

ax.legend()
plt.show()

# %%
# Compute rose21sd model without Ozone and plotting difference
O3AbsModel.model = ''
df_no_o3 = rte.execute()
df_no_o3 = df_no_o3.set_index(frq)
df['delta'] = df.tbtotal - df_no_o3.tbtotal

fig, ax = plt.subplots(1, 1, figsize=(12,8))
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('$\Delta {T_B}$ [K]')
df.delta.plot(ax=ax, figsize=(12,8), label='$\Delta {T_B}$ (rose21sd-rose21sd_w03)')
ax.legend()
plt.show()

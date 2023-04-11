"""
Performing Brightness Temperature calculation from ground with Ozone
====================================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.tb_spectrum.TbCloudRTE` method to calculate brightness temperature from ground (downwelling) with ozone.

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np

from pyrtlib.atmospheric_profiles import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.absorption_model import H2OAbsModel, O3AbsModel
from pyrtlib.utils import ppmv2gkg, mr2rh, ppmv_to_moleculesm3, import_lineshape, constants

atm = ['Tropical',
       'Midlatitude Summer',
       'Midlatitude Winter',
       'Subarctic Summer',
       'Subarctic Winter',
       'U.S. Standard']

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('${T_B}$ [K]')

z, p, d, t, md = atmp.gl_atm(atmp.US_STANDARD) # 'U.S. Standard'

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

rte = TbCloudRTE(z, p, t, rh, frq, ang, o3n)
rte.init_absmdl('rose20')
rte.satellite = False
H2OAbsModel.model = 'rose21sd'
H2OAbsModel.set_ll()
O3AbsModel.model = 'rose18'
O3AbsModel.set_ll()
df = rte.execute()

df = df.set_index(frq)
df.tbtotal.plot(ax=ax, linewidth=1, label='{} - {}'.format(atm[atmp.US_STANDARD], 'rose21sd'))

style = dict(size=20, color='gray', ha='center')
ax.text(22, 45, "${H_2O}$", **style)
ax.text(60, 255, "${O_2}}$", **style)
ax.text(119, 280, "${O_2}}$", **style)
ax.text(142, 100, "${O_3}}$", **style)
ax.text(183, 245, "${H_2O}$", **style)

def ghz_to_mm(ghz):
    f = ghz * 1e9
    c = constants('light')[0]
    return (c/f) * 1e3

def mm_to_ghz(mm):
    l = mm / 1e3
    c = constants('light')[0]
    return (c/l) / 1e9

secax = ax.secondary_xaxis('top', functions=(ghz_to_mm, mm_to_ghz))
secax.set_xlabel('$\lambda$ [mm]')

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

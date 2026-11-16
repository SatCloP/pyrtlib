"""
Performing Upwelling Brightness Temperature calculation using IGRA2 Upper Air Observations (with Extrapolation).
================================================================================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.tb_spectrum.TbCloudRTE` method to calculate brightness temperature from satellite (upwelling) using
# observations from IGRA2 Upper Air Archive and comparison of BT with the extrapoletd profile.

import numpy as np
from datetime import datetime

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})

from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.climatology import ProfileExtrapolation
from pyrtlib.utils import dewpoint2rh, to_kelvin
from pyrtlib.absorption_model import H2OAbsModel
from pyrtlib.apiwebservices import IGRAUpperAir

date = datetime(2020, 6, 1, 12)
station = 'SPM00008221'
df_igra2, header = IGRAUpperAir.request_data(date, station)

df_igra2 = df_igra2[df_igra2.pressure.notna() & 
                    df_igra2.temperature.notna() & 
                    df_igra2.dewpoint.notna() & 
                    df_igra2.height.notna()]

z, p, t = df_igra2.height.values / 1000, df_igra2.pressure.values, to_kelvin(df_igra2.temperature.values)

rh = dewpoint2rh(df_igra2.dewpoint, df_igra2.temperature).values

mdl = 'R21SD'
frq = np.arange(20, 201, 1)
nf = len(frq)

rte = TbCloudRTE(z, p, t, rh, frq)
rte.init_absmdl('R20')
H2OAbsModel.model = 'R21SD'
H2OAbsModel.set_ll()
df = rte.execute()
df = df.set_index(frq)

# %%
# Extrapolation of profile
ex = ProfileExtrapolation()
zz, pp, tt, rhh = ex.profile_extrapolation(header.latitude.values[0], 6, z, (p, t, rh))

rte = TbCloudRTE(zz, pp, tt, rhh, frq)
rte.init_absmdl('R20')
H2OAbsModel.model = 'R21SD'
H2OAbsModel.set_ll()
dff = rte.execute()
dff = dff.set_index(frq)

#%%
# Plotting
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
plt.suptitle("{}, {}, {} - {}".format(header.site_id.values[0], header.latitude.values[0], header.longitude.values[0], header.date.values[0]), y=0.96)
plt.title("IGRA2 UpperAir Radiosonde Archive", fontsize=10, ha='center')
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('${T_B}$ [K]')
df.tbtotal.plot(ax=ax, linewidth=2, label='{} - {}'.format(header.site_id.values[0], mdl))
dff.tbtotal.plot(ax=ax, linewidth=2, label='Extrap {} - {}'.format(header.site_id.values[0], mdl))
ax.grid(True, 'both')
ax.legend()
plt.show()

#%%
# Difference BT

df['delta'] = dff.tbtotal - df.tbtotal
df.delta.plot(linewidth=2, xlabel='Frequency [GHz]', ylabel='$\Delta T_B$ [K]', grid=True, figsize=(12, 8))

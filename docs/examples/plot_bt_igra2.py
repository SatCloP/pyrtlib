"""
Performing BT calculation from satellite using IGRA2 Upper Air Observations
=============================================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.main.BTCloudRTE` method to calculate brightness temperature from satellite using
# observations from Wyoming Upper Air Archive.

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np
from datetime import datetime

from pyrtlib.main import BTCloudRTE
from pyrtlib.utils import dewpoint2rh, import_lineshape
from pyrtlib.absmodel import H2OAbsModel
from pyrtlib.apiwebservices import IGRAUpperAir

date = datetime(2020, 6, 1, 12)
station = 'SPM00008221'
df_igra2, header = IGRAUpperAir.request_data(date, station)

df_igra2 = df_igra2[df_igra2.pressure.notna() & 
                    df_igra2.temperature.notna() & 
                    df_igra2.dewpoint.notna() & 
                    df_igra2.height.notna()]

z, p, t = df_igra2.height.values / 1000, df_igra2.pressure.values, df_igra2.temperature.values + 273.25

rh = dewpoint2rh(df_igra2.dewpoint, df_igra2.temperature).values

mdl = 'rose21sd'
ang = np.array([90.])
frq = np.arange(20, 201, 1)
nf = len(frq)

rte = BTCloudRTE(z, p, t, rh, frq, ang)
rte.init_absmdl('rose20')
H2OAbsModel.model = 'rose21sd'
H2OAbsModel.h2oll = import_lineshape('h2oll_{}'.format(H2OAbsModel.model))
df = rte.execute()
df = df.set_index(frq)

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
plt.suptitle("{}, {}, {} - {}".format(header.site_id.values[0], header.latitude.values[0], header.longitude.values[0], header.date.values[0]), y=0.96)
plt.title("IGRA2 UpperAir Radiosonde Archive", fontsize=10, ha='center')
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('${T_B}$ [K]')
df.tbtotal.plot(ax=ax, linewidth=2, label='{} - {}'.format(header.site_id.values[0], mdl))
ax.legend()
plt.show()

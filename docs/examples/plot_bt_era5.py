"""
Performing BT calculation form satellite using ERA5Reanalysis Observations
==========================================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.main.BTCloudRTE` method to calculate brightness temperature from satellite using
# observations from ERA5Reanalysis Reanalysis hourly pressure levels dataset.

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np
from pyrtlib.main import BTCloudRTE
from pyrtlib.utils import import_lineshape
from pyrtlib.absmodel import H2OAbsModel
from pyrtlib.apiwebservices import ERA5Reanalysis

# To request dataset from via ERA5Reanalysis API
# date = datetime(2020, 2, 22, 12)
# nc_file = ERA5Reanalysis.request_data(tempfile.gettempdir(), date, lonlat)

# Tito Scalo, Potenza, Italy
lonlat = (15.724447, 40.601019)
nc_file = 'tito_download_era5.nc'
z, p, t, rh, time = ERA5Reanalysis.read_data(nc_file, lonlat)

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
plt.title(
    "ERA5Reanalysis Reanalysis dataset (hourly pressure levels) {}".format(time[0].strftime(format='%Y-%m-%d %H:%M')),
    ha='center')
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('${T_B}$ [K]')
df.tbtotal.plot(ax=ax, linewidth=2, label='{} - {}'.format(lonlat, mdl))
ax.legend()
plt.show()

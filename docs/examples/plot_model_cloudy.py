"""
Performing Brightness Temperature calculation in cloudy condition
=================================================================
"""

# %%
# This example shows how to use the
# :py:meth:`pyrtlib.main.tb_cloud_rte` method to calculate brightness temperature from Satellite in cloudy condition

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FormatStrFormatter
plt.rcParams.update({'font.size': 15})
import numpy as np
np.seterr('raise')

from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import tb_cloud_rte
from pyrtlib.utils import ppmv2gkg, mr2rh

atm = ['Tropical',
       'Midlatitude Summer',
       'Midlatitude Winter',
       'Subarctic Summer',
       'Subarctic Winter',
       'U.S. Standard']

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

z, p, d, t, md = atmp.gl_atm(atmp.MIDLATITUDE_SUMMER)
gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
rh = mr2rh(p, t, gkg)[0] / 100

mdl = 'rose19sd'

ang = np.array([90.])
frq = np.arange(20, 201, 1)
nf = len(frq)

denliq = np.zeros(z.shape)
denice = np.zeros(z.shape)
cldh = np.empty((2, 2))

for i in [0, 1]:
    if i == 0:
        text_plot = 'clear-sky'
    else:
        # build a cloud
        ib = 1
        it = 3
        denliq[ib:it + 1] = 10 * np.ones((it - ib + 1))
        cldh[:, 0] = np.array([z[ib], z[it]])
        ib = 29
        it = 31
        denice[ib:it + 1] = 0.1 * np.ones((it - ib + 1))
        cldh[:, 1] = np.array([z[ib], z[it]])
        text_plot = 'cloudy'

    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('BT (K)')

    df = tb_cloud_rte(z, p, t, rh, denliq, denice, cldh, frq, ang,
                      absmdl=mdl,
                      ray_tracing=True,
                      from_sat=True,
                      cloudy=i)
    df = df.set_index(frq)
    df.tbtotal.plot(x=frq, ax=ax, linewidth=1, label='{} - {} ({})'.format(atm[atmp.MIDLATITUDE_SUMMER], mdl, text_plot))

ax.legend()
plt.show()

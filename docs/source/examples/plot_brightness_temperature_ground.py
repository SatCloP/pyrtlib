"""
Performing Brightness Temperature calculation from ground
=========================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.main.BTCloudRTE` method to calculate brightness temperature from ground

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np

from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh

atm = ['Tropical',
       'Midlatitude Summer',
       'Midlatitude Winter',
       'Subarctic Summer',
       'Subarctic Winter',
       'U.S. Standard']

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

for i in range(0, 6):
    z, p, d, t, md = atmp.gl_atm(i)
    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100

    mdl = 'rose19sd'

    ang = np.array([90.])
    frq = np.arange(20, 61, 1)
    nf = len(frq)

    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('BT (K)')

    rte = BTCloudRTE(z, p, t, rh, frq, ang)
    rte.satellite = False
    rte.init_absmdl(mdl)
    df = rte.execute()

    df = df.set_index(frq)
    df.tbtotal.plot(ax=ax, linewidth=1, label='{} - {}'.format(atm[i], mdl))

ax.legend()
plt.show()
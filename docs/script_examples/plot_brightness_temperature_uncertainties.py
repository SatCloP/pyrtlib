"""
Performing uncertainty on Brightness Temperature 
================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.main.BTCloudRTE` method to calculate uncertainty on brightness temperature
# with gamma_a as perturbation sptectroscopic parameter

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np

from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import BTCloudRTE
from pyrtlib.absmodel import H2OAbsModel, O3AbsModel
from pyrtlib.absmod_uncertainty import absmod_uncertainties_perturb
from pyrtlib.utils import ppmv2gkg, mr2rh

atm = ['Tropical',
       'Midlatitude Summer',
       'Midlatitude Winter',
       'Subarctic Summer',
       'Subarctic Winter',
       'U.S. Standard']

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.set_xlabel('Frequency [GHz]')
ax.set_ylabel('$\Delta {T_B}$ [K]')

for i in range(0, 6):

    z, p, d, t, md = atmp.gl_atm(i) # 'U.S. Standard'

    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100

    ang = np.array([90.])
    frq = np.arange(20, 151, 1)

    amu = absmod_uncertainties_perturb()

    rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
    rte.init_absmdl('uncertainty')
    df = rte.execute()
    df = df.set_index(frq)

    amu = absmod_uncertainties_perturb(['gamma_a'], 'max', index=1)

    rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
    rte.init_absmdl('uncertainty')
    df_gamma = rte.execute()
    df_gamma = df_gamma.set_index(frq)

    df['delta_max_gamma_a'] = df_gamma.tbtotal - df.tbtotal

    amu = absmod_uncertainties_perturb(['gamma_a'], 'min', index=1)

    rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
    rte.init_absmdl('uncertainty')
    df_gamma = rte.execute()
    df_gamma = df_gamma.set_index(frq)

    df['delta_min_gamma_a'] = df_gamma.tbtotal - df.tbtotal

    df.delta_max_gamma_a.plot(ax=ax, style='--', label='{}'.format(atm[i]))
    df.delta_min_gamma_a.plot(ax=ax, label='{}'.format(atm[i]))

# ax.legend()
plt.title("Perturbed parameters: $\gamma_a$")
plt.show()
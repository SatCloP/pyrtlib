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

z, p, d, t, md = atmp.gl_atm(5) # 'U.S. Standard'

gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
rh = mr2rh(p, t, gkg)[0] / 100

ang = np.array([90.])
frq = np.arange(120, 201, 1)

amu = absmod_uncertainties_perturb(['gamma_a', 'gamma_w'], 'non', index=2)

rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
rte.init_absmdl('uncertainty')
df = rte.execute()
df = df.set_index(frq)

amu = absmod_uncertainties_perturb(['gamma_a', 'gamma_w'], 'max', index=2)

rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
rte.init_absmdl('uncertainty')
df_gamma_a = rte.execute()
df_gamma_a = df_gamma_a.set_index(frq)

df['delta_max'] = df.tbtotal - df_gamma_a.tbtotal

amu = absmod_uncertainties_perturb(['gamma_a', 'gamma_w'], 'min', index=2)

rte = BTCloudRTE(z, p, t, rh, frq, ang, amu=amu)
rte.init_absmdl('uncertainty')
df_gamma_a = rte.execute()
df_gamma_a = df_gamma_a.set_index(frq)

df['delta_min'] = df.tbtotal - df_gamma_a.tbtotal

df.delta_max.plot(ax=ax, style='--', label='$\Delta {T_B}$ (uncertainty_max)')
df.delta_min.plot(ax=ax, label='$\Delta {T_B}$ (uncertainty_min)')

ax.legend()
plt.show()

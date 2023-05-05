"""
Performing uncertainty on Brightness Temperature 
================================================
"""

# %%
# This example shows how to use the
# :py:class:`pyrtlib.tb_spectrum.TbCloudRTE` method to calculate uncertainty on brightness temperature
# with gamma_a as perturbation sptectroscopic parameter

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
import numpy as np

from pyrtlib.atmospheric_profiles import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.absorption_model import H2OAbsModel, O3AbsModel
from pyrtlib.uncertainty import AbsModUncertainty
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

    amu = AbsModUncertainty.parameters_perturbation()

    rte = TbCloudRTE(z, p, t, rh, frq, ang, amu=amu)
    rte.init_absmdl('uncertainty')
    df = rte.execute()
    df = df.set_index(frq)

    amu = AbsModUncertainty.parameters_perturbation(['gamma_a'], 'max', index=1)

    rte = TbCloudRTE(z, p, t, rh, frq, ang, amu=amu)
    rte.init_absmdl('uncertainty')
    df_gamma = rte.execute()
    df_gamma = df_gamma.set_index(frq)

    df['delta_max_gamma_a'] = df_gamma.tbtotal - df.tbtotal

    amu = AbsModUncertainty.parameters_perturbation(['gamma_a'], 'min', index=1)

    rte = TbCloudRTE(z, p, t, rh, frq, ang, amu=amu)
    rte.init_absmdl('uncertainty')
    df_gamma = rte.execute()
    df_gamma = df_gamma.set_index(frq)

    df['delta_min_gamma_a'] = df_gamma.tbtotal - df.tbtotal

    df.delta_max_gamma_a.plot(ax=ax, style='--', label='{}'.format(atm[i]))
    df.delta_min_gamma_a.plot(ax=ax, label='{}'.format(atm[i]))

# ax.legend()
ax.grid(True, 'both')
plt.title("Perturbed parameters: $\gamma_a$")
plt.show()

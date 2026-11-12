"""
Computation of Weighting Functions
==================================
"""

# %%
# This example shows how to use the :py:class:`pyrtlib.weighting_functions.WeightingFunctions` method 
# to compute the weighting functions for the MWS channels for the U.S. standard atmospheric profile.

import numpy as np
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from pyrtlib.weighting_functions import WeightingFunctions
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.utils import ppmv2gkg, mr2rh, get_frequencies_sat

z, p, _, t, md = atmp.gl_atm(atmp.US_STANDARD)
gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
rh = mr2rh(p, t, gkg)[0] / 100

wf = WeightingFunctions(z, p, t, rh)
wf.frequencies = np.array([50.5, 53.2, 54.35, 54.9, 59.4, 58.825, 58.4])
wgt = wf.generate_wf()

wf.plot_wf(wgt, 'Downlooking', ylim=[0, 60], legend=True, figsize=(8, 6), dpi=100)

#%%
# As above but with the weighting functions computed in uplooking mode.

wf.satellite = False
wgt = wf.generate_wf()

wf.plot_wf(wgt, 'Uplooking', ylim=[0, 10], figsize=(8, 6), dpi=100)

#%%
# The weighting functions can also be computed for a different set of channels.
# The bandpass values are used to compute the weighting functions for the ATMS channels.
# The following code compute the weighting functions for the ATMS channels 5-15.

cf53 = 53.596
cf57 = 57.290344
frq = np.array([52.8, cf53-0.115, cf53+0.115, 54.4, 54.94, 55.5, cf57, 
       cf57-0.217, cf57+0.217, 
       cf57-0.3222-0.048, cf57-0.3222+0.048, cf57+0.3222-0.048, cf57+0.3222+0.048,
       cf57-0.3222-0.022, cf57-0.3222+0.022, cf57+0.3222-0.022, cf57+0.3222+0.022,
       cf57-0.3222-0.010, cf57-0.3222+0.010, cf57+0.3222-0.010, cf57+0.3222+0.010,
       cf57-0.3222-0.0045, cf57-0.3222+0.0045, cf57+0.3222-0.0045, cf57+0.3222+0.0045])

wf.satellite = True
wf.frequencies = frq
wf.bandpass = np.array([1, 2, 1, 1, 1, 1, 2, 4, 4, 4, 4])
wf.legend_labels = [f'Channel {i+5}' for i in range(len(wf.bandpass))]
wgt = wf.generate_wf()

wf.plot_wf(wgt, 'ATMS Channels 5-15', ylim=[0, 70], xlim=[0, 0.11], legend=True, figsize=(8, 6), dpi=100)

#%%
# The weighting functions can also be computed for a different set of frequencies.
# The following code compute the weighting functions for the MWS channels for a standard tropical atmosphere.
# for grouped frequencies.

wf.satellite = True
wf.frequencies = get_frequencies_sat('MWS')
wgt = wf.generate_wf()
wf.plot_wf_grouped(wgt, 'MWS Channels (grouped)', ylim=[0, 60], 
                  grouped_frequencies=[4, 9, 19, 1, 13],
                  grouped_labels=['23-52', '53-55', '57', '89', '164-229'], dpi=350)
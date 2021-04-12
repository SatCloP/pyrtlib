"""
This function provides the uncertainties affecting absorption model
coefficients I was able to find in litterature.
The baseline are the routines of Rosenkranz 2016 + modification to water-2-air by Koshelev et al 2015.

Hystory:
    * 2016/12/05 - Nico - First created
    * 2018/12/19 - Nico - Modified for adding 658 GHz line as in Rosenkranz 2018 (see if strcmp(mdl,'r18') at the end of WV parameters)
    * 2020/06/25 - Nico - Modified to account SD 183 parameters only (see if strcmp(mdl,'r20sd'))

.. note:: *Phil* The parameters that contribute most uncertainty are related to pressure-broadening: widths, mixing, & their T dependence. The intensities and their temperature dependence are obtained (at least in my model) from HITRAN2012. B2=(Ef+Ei)/2kTo,
    where Ef, Ei=upper and lower energy levels of the line; k=Boltzmann const; To=296K. Although HITRAN has uncertainties for six parameters, the energy level is not among them.  
    I expect that it does not contribute significant uncertainty to the model.

.. note:: *From Verdes et al 2005* For the pressure shifts, an estimated uncertainty for all lines belonging 
    to the same molecular species is assumed. This is a poor approximation as it is well known that the pressure shift is usually very different (with 
    possible change of sign) from one given transition to the next. Thus, the investigated pressure shift uncertainties were 300 kHz/Torr for 
    H2O and 50 kHz/Torr for O2 lines
        
.. note:: *From Makarov et al., 2011* "The fidelity of the new model to the spectrometer data is generally better than 2% between 54 and 65 GHz."
    "It is therefore important to distinguish between the uncertainty of the calculated absorption, and uncertainties in the values of the model?s coefficients. The coefficients are adjusted to fit measurements, and underestimation of some coefficients can be counterbalanced by excess values in other coefficients when absorption is calculated."

References
----------
.. [1] Koshelev et al., JQSRT, 112, 2704?2712, 2011
.. [2] Koshelev et al., JQSRT, 154, 24-27, 2015
.. [3] Koshelev et al., in preparation, 2018 (22 GHz)
.. [4] Koshelev et al., JQSRT, 196, 78?86, 2017 (118 GHz)
.. [5] Turner et al., TGRSS, 47, 10, 3326-37, 2009
.. [6] Tretyakov, JMS, 2016
"""

from collections import namedtuple
import numpy as np

UNKNOWN=0.0
MB2TORR=0.750062
USEKOSHELEV2017=1
USEKOSHELEV2017_WHAT='RAD'
C_CM=29979245800.0

FIELDS = (
    'value',
    'uncer',
    'units',
    'refer')

AMU = namedtuple('AMU', FIELDS)
AMU.__new__.__defaults__ = (None,) * len(AMU._fields)

ATTRIBUTES = {
    'w2a': AMU(
        value=1.2,
        uncer=0.05,
        units='unitless',
        refer='Koshelev et al., JQSRT, 2015'
    ),
    'con_CF': AMU(
        value=5.43e-10,
        uncer=UNKNOWN,
        units='1/(km*(mb^2*GHz^2))',
        refer='RTE Rosen98 routine'
    ),
    # The value used by P. Rosenkranz is slightly different 
    # The multipliers in Turner et al are referenced to Ros 1998 model. 
    # Most of the line parameters have been updated since then, causing small changes to absorption in the 
    # spectrum windows. Therefore, to keep the new model consistent with Turner's results, it's necessary to 
    # modify the multipliers slightly, but mainly in the foreign-continuum;
    'con_Cf_factr': AMU(
        # value=1.105,
        value=1.0976,
        uncer=np.sqrt(0.098 ** 2 + 0.03 ** 2),
        units='unitless',
        refer='Turner et al., TGRSS, 2009'
    ),
    'con_Cs': AMU(
        value=1.8e-08,
        uncer=UNKNOWN,
        units='1/(km*(mb^2*GHz^2))',
        refer='RTE Rosen98 routine'
    ),
    'con_Cs_factr': AMU(
        value=0.79,
        uncer=np.sqrt(0.17 ** 2 + 0.06 ** 2),
        units='unitless',
        refer='Turner et al., TGRSS, 2009'
    ),
    'con_Xf': AMU(
        value=3.0,
        uncer=0.8,
        units='unitless',
        refer='Tretyakov, JMS, 2016; Koshelev et al. 2011'
    ),
    'con_Xs': AMU(
        value=7.5,
        uncer=0.6,
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'FL': AMU(
        value=np.append([22.23507985,183.310087],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([5e-08,1e-06],np.tile(UNKNOWN,(1,13))),
        units='GHz',
        refer='Tretyakov, JMS, 2016'
    ),
    'S': AMU(
        value=np.append([4.39,774.6],np.tile(UNKNOWN,(1,13))) * 1e-25,
        uncer=np.append([0.043,7.7],np.tile(UNKNOWN,(1,13))) * 1e-25,
        units='cm/mol',
        refer='Tretyakov, JMS, 2016'
    ),
    # units in Ros. model are [Hz*cm^2]; the conversion factor is just speed of light in cm (P. Rosenkranz, personal communication)
    'S_ros': AMU(
        value=np.append([4.39,774.6],np.tile(UNKNOWN,(1,13))) * C_CM,
        uncer=np.append([0.043,7.7],np.tile(UNKNOWN,(1,13)))* C_CM,
        units='Hz*cm^2',
        refer='Tretyakov, JMS, 2016'
    ),
    'B2': AMU(
        value=np.array([2.144,0.668,6.179,1.541,1.048,3.595,5.048,1.405,3.597,2.379,2.852,0.159,2.391,0.396,1.441]),
        uncer=np.array([2.144,0.668,6.179,1.541,1.048,3.595,5.048,1.405,3.597,2.379,2.852,0.159,2.391,0.396,1.441]) / 100,
        units='unitless',
        refer='Rosenkranz, 2016 model + Tretyakov pers. comm.'
    ),
    'gamma_a': AMU(
        value=np.append([3.63,3.926],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.14,0.02],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016'
    ),
    'gamma_a_rad_k2017': AMU(
        value=np.append([3.598,3.926],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.022,0.02],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'gamma_a_video_k2017': AMU(
        value=np.append([3.397,3.926],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.079,0.02],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'gamma_a_combo_k2017': AMU(
        value=np.append([3.584,3.926],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.052,0.02],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    'gamma_w': AMU(
        value=np.append([17.6,19.7],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.5,0.5],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016'
    ),
    'gamma_w_rad_k2017': AMU(
        value=np.append([17.713,19.7],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.015,0.5],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'gamma_w_video_k2017': AMU(
        value=np.append([17.35,19.7],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.12,0.5],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'gamma_w_combo_k2017': AMU(
        value=np.append([17.707,19.7],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.045,0.5],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    'n_a': AMU(
        value=np.append([0.7,0.74],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.05,0.03],np.tile(UNKNOWN,(1,13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'n_w': AMU(
        value=np.append([1.2,0.78],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.5,0.08],np.tile(UNKNOWN,(1,13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'delta_a': AMU(
        value=np.append([0.0,-0.096],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.1,0.01],np.tile(UNKNOWN,(1,13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'delta_a_rad_k2017': AMU(
        value=np.append([-0.044,-0.096],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.005,0.01],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'delta_a_video_k2017': AMU(
        value=np.append([0.081,-0.096],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.015,0.01],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'delta_a_combo_k2017': AMU(
        value=np.append([-0.032,-0.096],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.038,0.01],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    'delta_w': AMU(
        value=np.append([-0.4,0.23],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([1.2,0.03],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016'
    ),
    'delta_w_rad_k2017': AMU(
        value=np.append([1.085,0.23],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.011,0.03],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'delta_w_video_k2017': AMU(
        value=np.append([1.53,0.23],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.06,0.03],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'delta_w_combo_k2017': AMU(
        value=np.append([1.099,0.23],np.tile(UNKNOWN,(1,13))),
        uncer=np.append([0.079,0.03],np.tile(UNKNOWN,(1,13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
}
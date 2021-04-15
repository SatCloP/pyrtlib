"""
This module provides the uncertainties affecting absorption model
coefficients I was able to find in litterature.
The baseline are the routines of Rosenkranz 2016 + modification to water-2-air by Koshelev et al 2015.

Example:

>>> from pyrtlib.absmod_uncertainty import AMU
>>> AMU['delta_a'].uncer
array([0.1 , 0.01, 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  , 0.  ,
       0.  , 0.  , 0.  , 0.  ])
>>> AMU['delta_a'].value
array([ 0.   , -0.096,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ])

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
import os
import numpy as np

from .utils import uncertainty_propagation

UNKNOWN = 0.0
MB2TORR = 0.750062
USEKOSHELEV2017 = 1
USEKOSHELEV2017_WHAT = 'RAD'
C_CM = 29979245800.0

PATH = os.path.dirname(os.path.abspath(__file__))
U = np.loadtxt(open(os.path.join(PATH, "lineshape", "u.csv"), "rb"), delimiter=",")
SIGMA = np.loadtxt(open(os.path.join(PATH, "lineshape", "sigma_widths_revised.csv"), "rb"), delimiter=",")

FIELDS = (
    'value',
    'uncer',
    'units',
    'refer')

AMU_NT = namedtuple('AMU_NT', FIELDS)
AMU_NT.__new__.__defaults__ = (None,) * len(AMU_NT._fields)

AMU = {
    'w2a': AMU_NT(
        value=1.2,
        uncer=0.05,
        units='unitless',
        refer='Koshelev et al., JQSRT, 2015'
    ),
    'con_CF': AMU_NT(
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
    'con_Cf_factr': AMU_NT(
        # value=1.105,
        value=1.0976,
        uncer=np.sqrt(0.098 ** 2 + 0.03 ** 2),
        units='unitless',
        refer='Turner et al., TGRSS, 2009'
    ),
    'con_Cs': AMU_NT(
        value=1.8e-08,
        uncer=UNKNOWN,
        units='1/(km*(mb^2*GHz^2))',
        refer='RTE Rosen98 routine'
    ),
    'con_Cs_factr': AMU_NT(
        value=0.79,
        uncer=np.sqrt(0.17 ** 2 + 0.06 ** 2),
        units='unitless',
        refer='Turner et al., TGRSS, 2009'
    ),
    'con_Xf': AMU_NT(
        value=3.0,
        uncer=0.8,
        units='unitless',
        refer='Tretyakov, JMS, 2016; Koshelev et al. 2011'
    ),
    'con_Xs': AMU_NT(
        value=7.5,
        uncer=0.6,
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'FL': AMU_NT(
        value=np.append([22.23507985, 183.310087], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([5e-08, 1e-06], np.tile(UNKNOWN, (1, 13))),
        units='GHz',
        refer='Tretyakov, JMS, 2016'
    ),
    'S': AMU_NT(
        value=np.append([4.39, 774.6], np.tile(UNKNOWN, (1, 13))) * 1e-25,
        uncer=np.append([0.043, 7.7], np.tile(UNKNOWN, (1, 13))) * 1e-25,
        units='cm/mol',
        refer='Tretyakov, JMS, 2016'
    ),
    # units in Ros. model are [Hz*cm^2]; the conversion factor is just speed of light in cm (P. Rosenkranz, personal communication)
    'S_ros': AMU_NT(
        value=np.append([4.39, 774.6], np.tile(UNKNOWN, (1, 13))) * C_CM,
        uncer=np.append([0.043, 7.7], np.tile(UNKNOWN, (1, 13))) * C_CM,
        units='Hz*cm^2',
        refer='Tretyakov, JMS, 2016'
    ),
    'B2': AMU_NT(
        value=np.array(
            [2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 0.159, 2.391, 0.396, 1.441]),
        uncer=np.array(
            [2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 0.159, 2.391, 0.396,
             1.441]) / 100,
        units='unitless',
        refer='Rosenkranz, 2016 model + Tretyakov pers. comm.'
    ),
    'gamma_a': AMU_NT(
        value=np.append([3.63, 3.926], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.14, 0.02], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016'
    ),
    'gamma_a_rad_k2017': AMU_NT(
        value=np.append([3.598, 3.926], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.022, 0.02], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'gamma_a_video_k2017': AMU_NT(
        value=np.append([3.397, 3.926], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.079, 0.02], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'gamma_a_combo_k2017': AMU_NT(
        value=np.append([3.584, 3.926], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.052, 0.02], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    'gamma_w': AMU_NT(
        value=np.append([17.6, 19.7], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.5, 0.5], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016'
    ),
    'gamma_w_rad_k2017': AMU_NT(
        value=np.append([17.713, 19.7], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.015, 0.5], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'gamma_w_video_k2017': AMU_NT(
        value=np.append([17.35, 19.7], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.12, 0.5], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'gamma_w_combo_k2017': AMU_NT(
        value=np.append([17.707, 19.7], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.045, 0.5], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    'n_a': AMU_NT(
        value=np.append([0.7, 0.74], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.05, 0.03], np.tile(UNKNOWN, (1, 13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'n_w': AMU_NT(
        value=np.append([1.2, 0.78], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.5, 0.08], np.tile(UNKNOWN, (1, 13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'delta_a': AMU_NT(
        value=np.append([0.0, -0.096], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.1, 0.01], np.tile(UNKNOWN, (1, 13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016'
    ),
    'delta_a_rad_k2017': AMU_NT(
        value=np.append([-0.044, -0.096], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.005, 0.01], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'delta_a_video_k2017': AMU_NT(
        value=np.append([0.081, -0.096], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.015, 0.01], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'delta_a_combo_k2017': AMU_NT(
        value=np.append([-0.032, -0.096], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.038, 0.01], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    'delta_w': AMU_NT(
        value=np.append([-0.4, 0.23], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([1.2, 0.03], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016'
    ),
    'delta_w_rad_k2017': AMU_NT(
        value=np.append([1.085, 0.23], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.011, 0.03], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
    ),
    'delta_w_video_k2017': AMU_NT(
        value=np.append([1.53, 0.23], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.06, 0.03], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
    ),
    'delta_w_combo_k2017': AMU_NT(
        value=np.append([1.099, 0.23], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.079, 0.03], np.tile(UNKNOWN, (1, 13))),
        units='MHz/Torr',
        refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
    ),
    # Shift to width ratio (this is used in Ros 2016 & 2017 models instead of delta_a/w)
    # In Ros 2016: SR = delta_a/gamma_a
    'SR': AMU_NT(
        value=np.append([0.0, 0.0], np.tile(UNKNOWN, (1, 13))),
        uncer=np.append([0.0, 0.0], np.tile(UNKNOWN, (1, 13))),
        units='unitless',
        refer='Tretyakov, JMS, 2016 - Eq. 3 and 5'
    ),
    'wv_nS': AMU_NT(
        value=2.5,
        uncer=2.5 * 0.005,
        units='unitless',
        refer='Rosenkranz, email of 2018/03/08, based on Gamache et al., 2017'
    ),
    # Line frequency FL (freq,GHz)
    # I take values from Ros 2016 and uncertainty from Tretyakov et al., 2005 (See Table 1) 
    # Note Table 1 shows the first 26 lines in Rosenk 2015, while 27+ is the the 28th of Rosenk 2015
    # For all higher freq lines I assume 17*1e-6 uncertainty (i.e. the largest value)
    'O2FL': AMU_NT(
        value=np.array(
            [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.591, 59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
             56.9682, 62.4112, 56.3634, 62.998, 55.7838, 63.5685, 55.2214, 64.1278, 54.6712, 64.6789, 54.13, 65.2241,
             53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, 52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.431,
             50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.763, 487.2493, 566.8956, 715.3929, 731.1866, 773.8395,
             834.1455, 895.071]),
        uncer=1e-06 * np.append(
            [7, 12, 6, 5, 5, 5, 4, 4, 4, 4, 4, 6, 4, 4, 5, 4, 6, 4, 5, 4, 7, 10, 15, 5, 17, 12, 12, 12],
            np.tile(17, (1, 21))),
        units='GHz',
        refer='Tretyakov et al., 2005'
    ),
    # From M. Tretyakov:
    # Line strengths were measured by Liebe, et al. 1977. 
    # We believe that HITRAN data are accurate to about 1#. 
    # We indirectly confirmed this by measuring N=1- intensity (paper in preparation). 
    # The coincidence with HITRAN is less than 0.5#. 
    'O2S': AMU_NT(
        value=np.array(
            [2.906e-15, 7.957e-16, 2.444e-15, 2.194e-15, 3.301e-15, 3.243e-15, 3.664e-15, 3.834e-15, 3.588e-15,
             3.947e-15, 3.179e-15, 3.661e-15, 2.59e-15, 3.111e-15, 1.954e-15, 2.443e-15, 1.373e-15, 1.784e-15,
             9.013e-16, 1.217e-15, 5.545e-16, 7.766e-16, 3.201e-16, 4.651e-16, 1.738e-16, 2.619e-16, 8.88e-17,
             1.387e-16, 4.272e-17, 6.923e-17, 1.939e-17, 3.255e-17, 8.301e-18, 1.445e-17, 3.356e-18, 6.049e-18,
             1.28e-18, 2.394e-18, 3.287e-17, 6.463e-16, 1.334e-17, 7.049e-15, 3.011e-15, 1.797e-17, 1.826e-15,
             2.193e-17, 1.153e-14, 3.974e-15, 2.512e-17]),
        uncer=np.array(
            [2.906e-15, 7.957e-16, 2.444e-15, 2.194e-15, 3.301e-15, 3.243e-15, 3.664e-15, 3.834e-15, 3.588e-15,
             3.947e-15, 3.179e-15, 3.661e-15, 2.59e-15, 3.111e-15, 1.954e-15, 2.443e-15, 1.373e-15, 1.784e-15,
             9.013e-16, 1.217e-15, 5.545e-16, 7.766e-16, 3.201e-16, 4.651e-16, 1.738e-16, 2.619e-16, 8.88e-17,
             1.387e-16, 4.272e-17, 6.923e-17, 1.939e-17, 3.255e-17, 8.301e-18, 1.445e-17, 3.356e-18, 6.049e-18,
             1.28e-18, 2.394e-18, 3.287e-17, 6.463e-16, 1.334e-17, 7.049e-15, 3.011e-15, 1.797e-17, 1.826e-15,
             2.193e-17, 1.153e-14, 3.974e-15, 2.512e-17]) / 100,
        units='cm2*Hz',
        refer='Tretyakov, Personal communication, 2016'
    ),
    'O2BE': AMU_NT(
        value=np.array(
            [0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621, 0.91, 0.91, 1.255, 1.255, 1.654,
             1.654, 2.109, 2.109, 2.618, 2.618, 3.182, 3.182, 3.8, 3.8, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819,
             6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651, 0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136,
             0.141, 0.145, 0.201]),
        uncer=np.array(
            [0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621, 0.91, 0.91, 1.255, 1.255, 1.654,
             1.654, 2.109, 2.109, 2.618, 2.618, 3.182, 3.182, 3.8, 3.8, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819,
             6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651, 0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136,
             0.141, 0.145, 0.201]) / 100,
        units='unitless',
        refer='Tretyakov, Personal communication, 2016'
    ),
    # TODO: UCM_O2 necessary
    # Line shape (i.e. pressure broadening parameters)
    # See Table 5 on Tretyakov et al., 2005 - a3 is (nearly) the same as W300 in Rosenk 2016 [GHz/(1e5 Pa)]
    # Only the first line differs, plus Ros has more lines at higher freq 
    # Table 5 does not report uncertainties; thus I use values on Table 1 and converts units
    # The uncertainties from table 1 have to be combined.
    # Note that Table 1 allows computation of combined uncertanity only for the first 20 (stronger) lines
    # In fact, only the first 20 lines have been investigated independently.
    # Since each line was measured independently, their uncertainties can be assumed as uncorrelated
    # Following Phil's suggestion, the combined uncertainty is obtained as sqrt(0.21*uO2^2 + 0.79*uN2^2).
    # These do not correspond to Phil's results listed in table 1 in his memo of July 11 (O2_off-diagonal.docx), because there they were computed as (0.21*uO2 + 0.79*uN2), which is not correct.
    # Units in Ros. model are [MHz/mb]; the conversion factor is 1/torr2mb, i.e. mb2torr
    # For other lines I here assume 0 uncertainty, as this is correlated and it's treated as O2gammaWL (weaker lines) below
    'O2BE': AMU_NT(
        value=np.array([0.01,0.014,0.083,0.083,0.207,0.207,0.387,0.387,0.621,0.621,0.91,0.91,1.255,1.255,1.654,1.654,2.109,2.109,2.618,2.618,3.182,3.182,3.8,3.8,4.474,4.474,5.201,5.201,5.983,5.983,6.819,6.819,7.709,7.709,8.653,8.653,9.651,9.651,0.019,0.048,0.045,0.044,0.049,0.084,0.145,0.136,0.141,0.145,0.201]),
        uncer= (0.25 * np.array([0.01,0.014,0.083,0.083,0.207,0.207,0.387,0.387,0.621,0.621,0.91,0.91,1.255,1.255,1.654,1.654,2.109,2.109,2.618,2.618,3.182,3.182,3.8,3.8,4.474,4.474,5.201,5.201,5.983,5.983,6.819,6.819,7.709,7.709,8.653,8.653,9.651,9.651,0.019,0.048,0.045,0.044,0.049,0.084,0.145,0.136,0.141,0.145,0.201])) / 100,
        units='unitless',
        refer='Tretyakov, Personal communication, 2016'
    ),
    'O2gamma': AMU_NT(
        value=np.array([1.688,1.703,1.513,1.491,1.415,1.408,1.353,1.339,1.295,1.292,1.262,1.263,1.223,1.217,1.189,1.174,1.134,1.134,1.089,1.088,1.037,1.038,0.996,0.996,0.955,0.955,0.906,0.906,0.858,0.858,0.811,0.811,0.764,0.764,0.717,0.717,0.669,0.669,2.78,1.64,1.64,1.64,1.6,1.6,1.6,1.6,1.62,1.47,1.47]),
        uncer= U[3:37],
        units='MHz/mb',
        refer='Tretyakov, JMS, 2005 + Rosenkranz Pers. Comm. 2017'
    ),
    # Line shape (i.e. pressure broadening parameters) for Weak Lines (WL)
    # From Phil (2017/09/22):
    # The width values for the weaker lines (21-38) in Tretyakov et al. (a3 in table 5) are extrapolated linearly
    # with quantum number N from the combined O2- and N2-broadening measured values of the stronger lines. 
    # That extrapolation introduces correlated uncertainties.
    # Per Phil suggestion (2017/07/17-18) I test 0.05 GHz/bar uncertainty for WL (21-38) and leave 0 uncertainty to the other
    # This is an arbitrary - likely overestimated - value: 0.05 GHz/bar = 0.05 MHz/mb
    # Phil recomputed the uncertainty considering linear regression (2017/10/02)
    # Per Phil suggestion (2017/10/03), after his new uncertainty calculations, these lines can be perturbed independently, so are treated above (O2gamma)
    # O2gamma_WL should only be used if one wants to compute the impact of uncertainties from all weaker lines together
    'O2gamma_WL': AMU_NT(
        value=np.array([1.688,1.703,1.513,1.491,1.415,1.408,1.353,1.339,1.295,1.292,1.262,1.263,1.223,1.217,1.189,1.174,1.134,1.134,1.089,1.088,1.037,1.038,0.996,0.996,0.955,0.955,0.906,0.906,0.858,0.858,0.811,0.811,0.764,0.764,0.717,0.717,0.669,0.669,2.78,1.64,1.64,1.64,1.6,1.6,1.6,1.6,1.62,1.47,1.47]),
        uncer= SIGMA[20:38, 2],
        units='MHz/mb',
        refer='Rosenkranz, pers. comm., 2017'
    ),
    # TODO: check better solution to referencing value within uncer
    'O2gamma_mmW': AMU_NT(
        value=np.array([1.688,1.703,1.513,1.491,1.415,1.408,1.353,1.339,1.295,1.292,1.262,1.263,1.223,1.217,1.189,1.174,1.134,1.134,1.089,1.088,1.037,1.038,0.996,0.996,0.955,0.955,0.906,0.906,0.858,0.858,0.811,0.811,0.764,0.764,0.717,0.717,0.669,0.669,2.78,1.64,1.64,1.64,1.6,1.6,1.6,1.6,1.62,1.47,1.47]),
        uncer= 0.0,
        units='MHz/mb',
        refer='Rosenkranz, pers. comm., 2017'
    ),
    # TODO: check better solution to referencing value within uncer
    'O2gamma_NL': AMU_NT(
        value=np.array([1.688,1.703,1.513,1.491,1.415,1.408,1.353,1.339,1.295,1.292,1.262,1.263,1.223,1.217,1.189,1.174,1.134,1.134,1.089,1.088,1.037,1.038,0.996,0.996,0.955,0.955,0.906,0.906,0.858,0.858,0.811,0.811,0.764,0.764,0.717,0.717,0.669,0.669,2.78,1.64,1.64,1.64,1.6,1.6,1.6,1.6,1.62,1.47,1.47]),
        uncer= 0.1 * SIGMA[20:38, 2],
        units='MHz/mb',
        refer='Rosenkranz, pers. comm., 2017'
    ),
    'Snr': AMU_NT(
        value=1.584e-17,
        uncer=1.584e-17 * 5 /100,
        units='Hz*cm2/GHz2',
        refer='Pickett et al., 1998, i.e. JPL line compilation'
    ),
    'WB300': AMU_NT(
        value=0.56,
        uncer=0.05,
        units='MHz/mb',
        refer='Rosenkranz, pers. comm., 2017'
    ),
    'X11': AMU_NT(
        value=0.785,
        uncer=0.035,
        units='unitless',
        refer='Makarov et al. JQSRT 2011 -> Makarov et al. JQSRT 2008'
    ),
    'X16': AMU_NT(
        value=0.765,
        uncer=0.011,
        units='unitless',
        refer='Koshelev et al., 2016'
    ),
    'X05': AMU_NT(
        value=0.8,
        uncer=0.05,
        units='unitless',
        refer='Tretyakov et al., 2005 + Koshelev et al., 2016'
    ),
    'APU': AMU_NT(
        value=1.0,
        uncer=1.0,
        units='unitless',
        refer='Makarov et al. JQSRT 2011'
    ),
    # TODO: check better solution to referencing value within uncer
    'Y300': AMU_NT(
        value=np.array([- 0.036,0.2547,- 0.3655,0.5495,- 0.5696,0.6181,- 0.4252,0.3517,- 0.1496,0.043,0.064,- 0.1605,0.2906,- 0.373,0.4169,- 0.4819,0.4963,- 0.5481,0.5512,- 0.5931,0.6212,- 0.6558,0.692,- 0.7208,0.7312,- 0.755,0.7555,- 0.7751,0.7914,- 0.8073,0.8307,- 0.8431,0.8676,- 0.8761,0.9046,- 0.9092,0.9416,- 0.9423,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
        uncer=1.0,
        units='1/bar == 1/1e5Pa =~ 1/atm',
        refer='Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017'
    ),
    # TODO: check better solution to referencing value within uncer
    'Y300_NL': AMU_NT(
        value=np.array([- 0.036,0.2547,- 0.3655,0.5495,- 0.5696,0.6181,- 0.4252,0.3517,- 0.1496,0.043,0.064,- 0.1605,0.2906,- 0.373,0.4169,- 0.4819,0.4963,- 0.5481,0.5512,- 0.5931,0.6212,- 0.6558,0.692,- 0.7208,0.7312,- 0.755,0.7555,- 0.7751,0.7914,- 0.8073,0.8307,- 0.8431,0.8676,- 0.8761,0.9046,- 0.9092,0.9416,- 0.9423,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
        uncer=1.0,
        units='1/bar == 1/1e5Pa =~ 1/atm',
        refer='Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017'
    ),
    'O2_V': AMU_NT(
        value=np.array([0.0079,- 0.0978,0.0844,- 0.1273,0.0699,- 0.0776,0.2309,- 0.2825,0.0436,- 0.0584,0.6056,- 0.6619,0.6451,- 0.6759,0.6547,- 0.6675,0.6135,- 0.6139,0.2952,- 0.2895,0.2654,- 0.259,0.375,- 0.368,0.5085,- 0.5002,0.6206,- 0.6091,0.6526,- 0.6393,0.664,- 0.6475,0.6729,- 0.6545,0.68,- 0.66,0.685,- 0.665,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
        uncer=1.0,
        units='1/bar == 1/1e5Pa =~ 1/atm',
        refer='Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017'
    ),
    'O2_V_NL': AMU_NT(
        value=np.array([0.0079,- 0.0978,0.0844,- 0.1273,0.0699,- 0.0776,0.2309,- 0.2825,0.0436,- 0.0584,0.6056,- 0.6619,0.6451,- 0.6759,0.6547,- 0.6675,0.6135,- 0.6139,0.2952,- 0.2895,0.2654,- 0.259,0.375,- 0.368,0.5085,- 0.5002,0.6206,- 0.6091,0.6526,- 0.6393,0.664,- 0.6475,0.6729,- 0.6545,0.68,- 0.66,0.685,- 0.665,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]),
        uncer=1.0,
        units='1/mb',
        refer='Rosenkranz, pers. comm., 2017'
    ),
    'O2_nS': AMU_NT(
        value=2.0,
        uncer=1.0 * 0.001,
        units='unitless',
        refer='Rosenkranz, email of 2017/10/06, based on Gamache et al., 2017'
    ),
}

np_zeros = np.zeros(AMU['O2gamma'].value.shape)
np_zeros[38:49] = AMU['O2gamma'].value[38:49] * 0.1
AMU['O2gamma_mmW'] = AMU['O2gamma_mmW']._replace(uncer=np_zeros)
np_zeros = np.zeros(AMU['O2gamma'].value.shape)
np_zeros[34:49] = AMU['O2gamma'].value[34:49] * 0.1
AMU['O2gamma_NL'] = AMU['O2gamma_NL']._replace(uncer=np_zeros)
np_zeros = np.zeros(AMU['Y300'].value.shape)
np_zeros[0:34] = U[37:71]
AMU['Y300'] = AMU['Y300']._replace(uncer=np_zeros)
np_zeros = np.zeros(AMU['Y300'].value.shape)
np_zeros[34:49] = AMU['Y300'].value[34:49] * 0.2
AMU['Y300_NL'] = AMU['Y300_NL']._replace(uncer=np_zeros)
np_zeros = np.zeros(AMU['O2_V'].value.shape)
np_zeros[0:34] = U[71:105]
AMU['O2_V'] = AMU['O2_V']._replace(uncer=np_zeros)
np_zeros = np.zeros(AMU['O2_V'].value.shape)
np_zeros[34:49] = AMU['O2_V'].value[34:49] * 0.2
AMU['O2_V_NL'] = AMU['O2_V_NL']._replace(uncer=np_zeros)

for i in range(0, 2):
    AMU['SR'].value[i], \
    AMU['SR'].uncer[i] = uncertainty_propagation('A/B',
                                                 AMU['delta_a'].value[i],
                                                 AMU['gamma_a'].value[i],
                                                 AMU['delta_a'].uncer[i],
                                                 AMU['gamma_a'].uncer[i])

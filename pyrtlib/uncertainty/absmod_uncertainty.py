"""
.. note:: *Phil* The parameters that contribute most uncertainty are related to pressure-broadening: widths, mixing, & their T dependence. The intensities and their temperature dependence are obtained (at least in my model) from HITRAN2012. B2=(Ef+Ei)/2kTo,
    where Ef, Ei=upper and lower energy levels of the line; k=Boltzmann const; To=296K. Although HITRAN has uncertainties for six parameters, the energy level is not among them.  
    I expect that it does not contribute significant uncertainty to the model.

.. note:: *From Verdes et al 2005* For the pressure shifts, an estimated uncertainty for all lines belonging 
    to the same molecular species is assumed. This is a poor approximation as it is well known that the pressure shift is usually very different (with 
    possible change of sign) from one given transition to the next. Thus, the investigated pressure shift uncertainties were 300 kHz/Torr for 
    H2O and 50 kHz/Torr for O2 lines
        
.. note:: *From Makarov et al., 2011* "The fidelity of the new model to the spectrometer data is generally better than 2% between 54 and 65 GHz."
    "It is therefore important to distinguish between the uncertainty of the calculated absorption, and uncertainties in the values of the model?s coefficients. The coefficients are adjusted to fit measurements, and underestimation of some coefficients can be counterbalanced by excess values in other coefficients when absorption is calculated."

"""

import os
from copy import deepcopy
from dataclasses import dataclass
from typing import Tuple, Optional, Dict
from dataclasses import dataclass, field

import numpy as np
# import scipy.interpolate as si
from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, O3AbsModel
from pyrtlib.utils import constants

# TVEC = np.loadtxt(
#     open(os.path.join(PATH, "tbd", "Tvec.csv"), "rb"), delimiter=",")
# PVEC = np.loadtxt(
#     open(os.path.join(PATH, "tbd", "Pvec.csv"), "rb"), delimiter=",")
# FVEC = np.loadtxt(
#     open(os.path.join(PATH, "tbd", "Fvec.csv"), "rb"), delimiter=",")
# PR_EXT = np.loadtxt(
#     open(os.path.join(PATH, "tbd", "PR_ext.csv"), "rb"), delimiter=",")


@dataclass
class SpectroscopicParameter:
    """Absorption model uncertainties for the spectroscopic parameters

    Example:
        .. code-block:: python

            >>> from pyrtlib.uncertainty import SpectroscopicParameter
            >>> parameters = SpectroscopicParameter.oxygen_parameters('R18')
            >>> parameters['O2gamma_WL'].value
            array([1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295,
                1.292, 1.262, 1.263, 1.223, 1.217, 1.189, 1.174, 1.134, 1.134,
                1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906,
                0.906, 0.858, 0.858, 0.811, 0.811, 0.764, 0.764, 0.717, 0.717,
                0.669, 0.669, 2.78 , 1.64 , 1.64 , 1.64 , 1.6  , 1.6  , 1.6  ,
                1.6  , 1.62 , 1.47 , 1.47 ])
            >>> parameters['O2gamma_WL'].uncer
            array([0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
                0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
                0.   , 0.   , 0.012, 0.012, 0.015, 0.015, 0.017, 0.017, 0.019,
                0.019, 0.021, 0.021, 0.024, 0.024, 0.026, 0.026, 0.028, 0.028,
                0.031, 0.031, 0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
                0.   , 0.   , 0.   , 0.   ])
            >>> parameters['O2gamma_WL'].units
            'MHz/mb'
            >>> parameters['O2gamma_WL'].refer
            'Rosenkranz, pers. comm., 2017'
    """
    value: np.ndarray
    """The value associated to the parameter
    """
    uncer: np.ndarray
    """The uncertainty of the parameter
    """
    units: Optional[str] = ''
    """The units of the parameter
    """
    refer: Optional[str] = ''
    """The reference of the parameter
    """
    name: Optional[str] = ''
    """The name or description of the parameter
    """

    def _initialize() -> Dict:
        SpectroscopicParameter._parameters()

    @staticmethod
    def water_parameters(model: str) -> Dict:
        """This method is used for uncertainty analysis and returns the dictionary
        with the whole spectroscopic parameters for :math:`H_2O`.

        Args:
            model (str): The absorption model.

        Returns:
            Dict: The :math:`H_2O` spectroscopic parameters dictionary
        """
        H2OAbsModel.model = model
        H2OAbsModel.set_ll()
        ll = H2OAbsModel.h2oll

        h2o_sp = {
            "con_Cf": SpectroscopicParameter(value=ll.cf, uncer=.0, units = '1/(km*(mb^2*GHz^2))', name='Foreign induced broadening coefficient'),
            "con_Cs": SpectroscopicParameter(value=ll.cs, uncer=.0, units = '1/(km*(mb^2*GHz^2))', name='Self induced broadening coefficient'),
            "con_Cf_factr": SpectroscopicParameter(value=1.11, uncer=.0, name=''),
            "con_Cs_factr": SpectroscopicParameter(value=0.79, uncer=.0, name=''),
            "con_Xf": SpectroscopicParameter(value=ll.xcf, uncer=.0, units = 'unitless', name='Foreign broadening temperature dependence exponents'),
            "con_Xs": SpectroscopicParameter(value=ll.xcs, uncer=.0, units = 'unitless', name='Self broadening temperature dependence exponents'),

            "FL": SpectroscopicParameter(value=ll.fl, uncer=.0, units = 'GHz', name='Central frequency line'),
            "S": SpectroscopicParameter(value=ll.s1, uncer=np.zeros(len(ll.fl)), units = 'Hz*cm^2', name='Line intensity (or strength)'),
            "B2": SpectroscopicParameter(value=ll.b2, uncer=.0, units = 'unitless', name='Temperature coefficient of intensity'),

            "gamma_a": SpectroscopicParameter(value=ll.w0*1e3, units='MHz/mb', uncer=np.zeros(len(ll.fl)), name='Air induced broadening coefficients'),
            "gamma_w": SpectroscopicParameter(value=ll.w0s*1e3, units = 'MHz/mb', uncer=np.zeros(len(ll.fl)), name='Water induced broadening coefficients'),
            "n_a": SpectroscopicParameter(value=ll.x, uncer=.0, units = 'unitless', name='Temperature-dependence exponent of air'),
            "n_w": SpectroscopicParameter(value=ll.xs, uncer=.0, units = 'unitless', name='Temperature-dependence exponent of water'),

            "wv_nS": SpectroscopicParameter(value=.0, uncer=.0, name='Intensity temperature-dependence exponent'),
        }
        if H2OAbsModel.model > 'R17':
            h2o_sp_ = {
                "n_da": SpectroscopicParameter(value=ll.xh, uncer=.0, name='T-exponent of air-shifting'),
                "n_dw": SpectroscopicParameter(value=ll.xhs, uncer=.0, name='T-exponent of water-shifting'),
                "delta_a": SpectroscopicParameter(value=ll.sh*1e3, uncer=np.zeros(len(ll.fl)), units = 'MHz/mb', name='Air induced shifting coefficient'),
                "delta_w": SpectroscopicParameter(value=ll.shs*1e3, uncer=np.zeros(len(ll.fl)), units = 'MHz/mb', name='Water induced shifting coefficient'),
            }
        else:
            h2o_sp_ = {"SR": SpectroscopicParameter(value=ll.sr, uncer=np.zeros(len(ll.fl)), units='unitless', name='Shift to width ratio'),}

        h2o_sp = {**h2o_sp, **h2o_sp_}

        return h2o_sp

    @staticmethod
    def oxygen_parameters(model: str) -> Dict:
        """This method is used for uncertainty analysis and returns the dictionary
        with the whole spectroscopic parameters for :math:`O_2`.

        Args:
            model (str): The absorption model.

        Returns:
            Dict: The :math:`O_2` spectroscopic parameters dictionary
        """
        O2AbsModel.model = model
        O2AbsModel.set_ll()
        ll = O2AbsModel.o2ll

        o2_sp = {
            "O2FL": SpectroscopicParameter(value=ll.f, uncer=.0, units = 'GHz', name='Line frequency'),
            "O2S": SpectroscopicParameter(value=ll.s300, uncer=.0, units = 'cm2*Hz', name='Line intensity (or strength)'),
            "O2BE": SpectroscopicParameter(value=ll.be, uncer=.0, units = 'unitless', name='Temperature exponent for intensity'),
            "O2gamma": SpectroscopicParameter(value=ll.w300, uncer=np.zeros(len(ll.f)), units = 'MHz/mb', name='Self broadening temperature dependence exponents'),
            "O2gamma_NL": SpectroscopicParameter(value=ll.w300, uncer=.0, units='MHz/mb', name='Self broadening temperature dependence exponents for neglected lines (NL)'),
            "Snr": SpectroscopicParameter(value=1.6e-17, uncer=.0, units = 'Hz*cm2/GHz2', name='Line intensity of non-resonant pseudo-line'),
            "WB300": SpectroscopicParameter(value=ll.wb300, uncer=.0, units = 'MHz/mb', name='Pressure broadening of non-resonant pseudo-line'),
            "X11": SpectroscopicParameter(value=0.785, uncer=.0, units = 'unitless', name='Temperature dependence of broadening coefficient'),
            "X16": SpectroscopicParameter(value=0.765, uncer=.0, units = 'unitless', name='Temperature dependence of broadening coefficient'),
            "X05": SpectroscopicParameter(value=0.80, uncer=.0, units = 'unitless', name='Temperature dependence of broadening coefficient'),
            "Y300": SpectroscopicParameter(value=ll.y300, uncer=np.zeros(len(ll.f)), units = '1/bar == 1/1e5Pa =~ 1/atm', name='Line mixing coefficients (single lines)'),
            "Y300_NL": SpectroscopicParameter(value=ll.y300, uncer=.0, units = '1/mb', name='Line mixing coefficients for neglected lines (NL)'),
            "O2_V": SpectroscopicParameter(value=ll.v, uncer=np.zeros(len(ll.f)), units = '1/bar == 1/1e5Pa =~ 1/atm', name='O2 line mixing temperature dependence'),
            "O2_V_NL": SpectroscopicParameter(value=ll.v, uncer=.0, units = '1/mb', name='O2 line mixing temperature dependence for neglected lines (NL)'),
            "O2_nS": SpectroscopicParameter(value=2.0, uncer=.0, units = 'unitless',  name='Intensity temperature-dependence exponent'),
            "w2a": SpectroscopicParameter(value=1.2, uncer=.0, units = 'unitless', name='Water-to-air broadening ratio'),
        }
        if O2AbsModel.model > 'R19':
            o2_sp_ = {
                "y0": SpectroscopicParameter(value=ll.y0, uncer=.0, name=''),
                "y1": SpectroscopicParameter(value=ll.y1, uncer=.0, name=''),
                "dnu0": SpectroscopicParameter(value=ll.dnu0, uncer=.0, name=''),
                "dnu1": SpectroscopicParameter(value=ll.dnu1, uncer=.0, name=''),
            }
        else:
            o2_sp_ = {}

        o2_sp = {**o2_sp, **o2_sp_}

        return o2_sp
    
    @staticmethod
    def ozono_parameters(model: str) -> Dict:
        """This method is used for uncertainty analysis and returns the dictionary
        with the whole spectroscopic parameters for :math:`O_3`.

        Args:
            model (str): The absorption model.

        Returns:
            Dict: The Ozono (:math:`O_3`) spectroscopic parameters dictionary
        """
        O3AbsModel.model = model
        O3AbsModel.set_ll()
        ll = O3AbsModel.o3ll

        o3_sp = {
            "O3_FL": SpectroscopicParameter(value=ll.fl, uncer=.0, units='GHz', name='Line frequency'),
            "O3_S1": SpectroscopicParameter(value=ll.s1, uncer=.0, units = 'Hz*cm^2', name='Line intensity (or strength)'),
            "O3_B": SpectroscopicParameter(value=ll.b, uncer=.0, units = 'unitless', name='Temperature Coefficient of intensity'),
            "O3_W": SpectroscopicParameter(value=ll.w, uncer=.0, units = 'GHz/mb', name='Air-pressure broadening'),
            "O3_X": SpectroscopicParameter(value=ll.x, uncer=.0, units = 'unitless', name='Temperature exponent of air broadening'),
            "O3_SR": SpectroscopicParameter(value=ll.sr, uncer=.0, units = 'unitless', name='Shift-to-width ratio'),
        }

        return o3_sp

    @staticmethod
    def set_parameters(SP: dict) -> None:
        """Set new values and uncertainties to spectroscopic parameters.

        Args:
            SP (dict): The new dictionary with the spectroscopic parameters

        Example:
            >>> parameters = SpectroscopicParameter.water_parameters('R17')
            >>> parameters['gamma_a'].value[0] = 2.688
            >>> parameters['gamma_a'].uncer[0] = 0.039
            >>> SpectroscopicParameter.set_parameters(parameters)

        """
        SpectroscopicParameter.SP = SP

    @staticmethod
    def _parameters() -> Dict:
        """Returns a mutable dictionary of all the spectroscopic parameters.

        Returns:
            dict: The spectroscopic parameters dictionary

        Example:
            >>> from pyrtlib.uncertainty import SpectroscopicParameter
            >>> parameters = SpectroscopicParameter.parameters()
            >>> parameters['w2a'].value
            1.2

        New value may be added to parameters using :py:func:`~pyrtlib.uncertainty.SpectroscopicParameter` class as following

        Example:
            >>> parameters['con_Xs'] = SpectroscopicParameter(2.3, 0.001, 'unitless', 'Tretyakov, JMS, 2016')
            >>> parameters['con_Xs'].value
            2.3

        Also, existing parameters may be modified as following

        Example:
            >>> parameters['w2a'].value = 1.333
            >>> parameters['w2a'].value
            1.333

        Note:
            If new value will be added or modified it is necessary to save the new values by calling
            the :py:func:`~pyrtlib.uncertainty.SpectroscopicParameter.set_parameters` function.
        """
        UNKNOWN = 0.0
        FL_LENGHT = 13
        MB2TORR = 0.750062
        USEKOSHELEV2017 = True
        USEKOSHELEV2017_WHAT = 'rad'
        C_CM = constants('light')[0] * 100

        PATH = os.path.dirname(os.path.abspath(__file__))
        U = np.loadtxt(
            open(os.path.join(PATH, "tbd", "u.csv"), "rb"), delimiter=",")
        SIGMA = np.loadtxt(open(os.path.join(
            PATH, "tbd", "sigma_widths_revised.csv"), "rb"), delimiter=",")

        SP = {
            'w2a': SpectroscopicParameter(
                value=1.2,
                uncer=0.05,
                units='unitless',
                refer='Koshelev et al., JQSRT, 2015'
            ),
            'con_Cf': SpectroscopicParameter(
                value=5.96e-10,
                uncer=5.5e-11,
                units='1/(km*(mb^2*GHz^2))',
                refer='RTE Rosen98 routine'
            ),
            # The value used by P. Rosenkranz is slightly different
            # The multipliers in Turner et al are referenced to Ros 1998 model.
            # Most of the line parameters have been updated since then, causing small changes to absorption in the
            # spectrum windows. Therefore, to keep the new model consistent with Turner's results, it's necessary to
            # modify the multipliers slightly, but mainly in the foreign-continuum;
            'con_Cf_factr': SpectroscopicParameter(
                value=1.11,
                uncer=np.sqrt(0.098 ** 2 + 0.03 ** 2),
                units='unitless',
                refer='Turner et al., TGRSS, 2009'
            ),
            'con_Cs': SpectroscopicParameter(
                value=1.42e-8,
                uncer=3.2e-9,
                units='1/(km*(mb^2*GHz^2))',
                refer='RTE Rosen98 routine'
            ),
            'con_Cs_factr': SpectroscopicParameter(
                value=0.79,
                uncer=np.sqrt(0.17 ** 2 + 0.06 ** 2),
                units='unitless',
                refer='Turner et al., TGRSS, 2009'
            ),
            'con_Xf': SpectroscopicParameter(
                value=3.0,
                uncer=0.8,
                units='unitless',
                refer='Tretyakov, JMS, 2016; Koshelev et al. 2011'
            ),
            'con_Xs': SpectroscopicParameter(
                value=7.5,
                uncer=0.6,
                units='unitless',
                refer='Tretyakov, JMS, 2016'
            ),
            'FL': SpectroscopicParameter(
                value=np.append([22.23507985, 183.310087],
                                np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append(
                    [5e-08, 1e-06], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='GHz',
                refer='Tretyakov, JMS, 2016'
            ),
            'S': SpectroscopicParameter(
                value=np.append([4.39, 774.6], np.tile(
                    UNKNOWN, (1, FL_LENGHT))) * 1e-25,
                uncer=np.append([0.043, 7.7], np.tile(
                    UNKNOWN, (1, FL_LENGHT))) * 1e-25,
                units='cm/mol',
                refer='Tretyakov, JMS, 2016'
            ),
            # units in Ros. model are [Hz*cm^2]; the conversion factor is just speed of light in cm (P. Rosenkranz, personal communication)
            'S_ros': SpectroscopicParameter(
                value=np.append([4.39, 774.6], np.tile(
                    UNKNOWN, (1, FL_LENGHT))) * C_CM,
                uncer=np.append([0.043, 7.7], np.tile(
                    UNKNOWN, (1, FL_LENGHT))) * C_CM,
                units='Hz*cm^2',
                refer='Tretyakov, JMS, 2016'
            ),
            'B2': SpectroscopicParameter(
                value=np.array(
                    [2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 0.159, 2.391, 0.396, 1.441]),
                uncer=np.array(
                    [2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 0.159, 2.391, 0.396,
                     1.441]) / 100,
                units='unitless',
                refer='Rosenkranz, 2016 model + Tretyakov pers. comm.'
            ),
            'gamma_a': SpectroscopicParameter(
                value=np.append([3.63, 3.926], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.14, 0.02], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016'
            ),
            'gamma_a_rad_k2017': SpectroscopicParameter(
                value=np.append([3.598, 3.926], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.022, 0.02], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
            ),
            'gamma_a_video_k2017': SpectroscopicParameter(
                value=np.append([3.397, 3.926], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.079, 0.02], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
            ),
            'gamma_a_combo_k2017': SpectroscopicParameter(
                value=np.append([3.584, 3.926], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.052, 0.02], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
            ),
            'gamma_w': SpectroscopicParameter(
                value=np.append([17.6, 19.7], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.5, 0.5], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016'
            ),
            'gamma_w_rad_k2017': SpectroscopicParameter(
                value=np.append([17.713, 19.7], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.015, 0.5], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
            ),
            'gamma_w_video_k2017': SpectroscopicParameter(
                value=np.append([17.35, 19.7], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.12, 0.5], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
            ),
            'gamma_w_combo_k2017': SpectroscopicParameter(
                value=np.append([17.707, 19.7], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.045, 0.5], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
            ),
            'n_a': SpectroscopicParameter(
                value=np.append([0.7, 0.74], np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.05, 0.03], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='unitless',
                refer='Tretyakov, JMS, 2016'
            ),
            'n_w': SpectroscopicParameter(
                value=np.append([1.2, 0.78], np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.5, 0.08], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='unitless',
                refer='Tretyakov, JMS, 2016'
            ),
            'delta_a': SpectroscopicParameter(
                value=np.append(
                    [0.0, -0.096], np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.1, 0.01], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='unitless',
                refer='Tretyakov, JMS, 2016'
            ),
            'delta_a_rad_k2017': SpectroscopicParameter(
                value=np.append([-0.044, -0.096],
                                np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.005, 0.01], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
            ),
            'delta_a_video_k2017': SpectroscopicParameter(
                value=np.append([0.081, -0.096],
                                np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.015, 0.01], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
            ),
            'delta_a_combo_k2017': SpectroscopicParameter(
                value=np.append([-0.032, -0.096],
                                np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.038, 0.01], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
            ),
            'delta_w': SpectroscopicParameter(
                value=np.append(
                    [-0.4, 0.23], np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([1.2, 0.03], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016'
            ),
            'delta_w_rad_k2017': SpectroscopicParameter(
                value=np.append([1.085, 0.23], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.011, 0.03], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, RAD)'
            ),
            'delta_w_video_k2017': SpectroscopicParameter(
                value=np.append([1.53, 0.23], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.06, 0.03], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Video)'
            ),
            'delta_w_combo_k2017': SpectroscopicParameter(
                value=np.append([1.099, 0.23], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.079, 0.03], np.tile(
                    UNKNOWN, (1, FL_LENGHT))),
                units='MHz/Torr',
                refer='Tretyakov, JMS, 2016 + Koshelev et al. 2018 (22 GHz, Combined)'
            ),
            # Shift to width ratio (this is used in Ros 2016 & 2017 models instead of delta_a/w)
            # In Ros 2016: SR = delta_a/gamma_a
            'SR': SpectroscopicParameter(
                value=np.append([0.0, 0.0], np.tile(UNKNOWN, (1, FL_LENGHT))),
                uncer=np.append([0.0, 0.0], np.tile(UNKNOWN, (1, FL_LENGHT))),
                units='unitless',
                refer='Tretyakov, JMS, 2016 - Eq. 3 and 5'
            ),
            'wv_nS': SpectroscopicParameter(
                value=2.5,
                uncer=2.5 * 0.005,
                units='unitless',
                refer='Rosenkranz, email of 2018/03/08, based on Gamache et al., 2017'
            ),
            # Line frequency FL (freq,GHz)
            # I take values from Ros 2016 and uncertainty from Tretyakov et al., 2005 (See Table 1)
            # Note Table 1 shows the first 26 lines in Rosenk 2015, while 27+ is the the 28th of Rosenk 2015
            # For all higher freq lines I assume 17*1e-6 uncertainty (i.e. the largest value)
            'O2FL': SpectroscopicParameter(
                value=np.array(
                    [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.591, 59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
                     56.9682, 62.4112, 56.3634, 62.998, 55.7838, 63.5685, 55.2214, 64.1278, 54.6712, 64.6789, 54.13, 65.2241,
                     53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, 52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.431,
                     50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.763, 487.2493, 566.8956, 715.3929, 731.1866, 773.8395,
                     834.1455, 895.071]),
                uncer=1e-06 * np.append(
                    [7, 12, 6, 5, 5, 5, 4, 4, 4, 4, 4, 6, 4, 4, 5,
                     4, 6, 4, 5, 4, 7, 10, 15, 5, 17, 12, 12, 12],
                    np.tile(17, (1, 21))),
                units='GHz',
                refer='Tretyakov et al., 2005'
            ),
            # From M. Tretyakov:
            # Line strengths were measured by Liebe, et al. 1977.
            # We believe that HITRAN data are accurate to about 1#.
            # We indirectly confirmed this by measuring N=1- intensity (paper in preparation).
            # The coincidence with HITRAN is less than 0.5#.
            'O2S': SpectroscopicParameter(
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
            'O2BE': SpectroscopicParameter(
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
            'O2BE': SpectroscopicParameter(
                value=np.array(
                    [0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621, 0.91, 0.91, 1.255, 1.255, 1.654,
                     1.654, 2.109, 2.109, 2.618, 2.618, 3.182, 3.182, 3.8, 3.8, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819,
                     6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651, 0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136,
                     0.141, 0.145, 0.201]),
                uncer=(0.25 * np.array(
                    [0.01, 0.014, 0.083, 0.083, 0.207, 0.207, 0.387, 0.387, 0.621, 0.621, 0.91, 0.91, 1.255, 1.255, 1.654,
                     1.654, 2.109, 2.109, 2.618, 2.618, 3.182, 3.182, 3.8, 3.8, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819,
                     6.819, 7.709, 7.709, 8.653, 8.653, 9.651, 9.651, 0.019, 0.048, 0.045, 0.044, 0.049, 0.084, 0.145, 0.136,
                     0.141, 0.145, 0.201])) / 100,
                units='unitless',
                refer='Tretyakov, Personal communication, 2016'
            ),
            'O2gamma': SpectroscopicParameter(
                value=np.array(
                    [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189,
                     1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858,
                     0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 2.78, 1.64, 1.64, 1.64, 1.6, 1.6, 1.6, 1.6, 1.62,
                     1.47, 1.47]),
                uncer=U[3:37],
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
            'O2gamma_WL': SpectroscopicParameter(
                value=np.array(
                    [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189,
                     1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858,
                     0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 2.78, 1.64, 1.64, 1.64, 1.6, 1.6, 1.6, 1.6, 1.62,
                     1.47, 1.47]),
                uncer=SIGMA[20:38, 2],
                units='MHz/mb',
                refer='Rosenkranz, pers. comm., 2017'
            ),
            # TODO: check better solution to referencing value within uncer
            'O2gamma_mmW': SpectroscopicParameter(
                value=np.array(
                    [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189,
                     1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858,
                     0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 2.78, 1.64, 1.64, 1.64, 1.6, 1.6, 1.6, 1.6, 1.62,
                     1.47, 1.47]),
                uncer=0.0,
                units='MHz/mb',
                refer='Rosenkranz, pers. comm., 2017'
            ),
            # TODO: check better solution to referencing value within uncer
            'O2gamma_NL': SpectroscopicParameter(
                value=np.array(
                    [1.688, 1.703, 1.513, 1.491, 1.415, 1.408, 1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217, 1.189,
                     1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038, 0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858,
                     0.811, 0.811, 0.764, 0.764, 0.717, 0.717, 0.669, 0.669, 2.78, 1.64, 1.64, 1.64, 1.6, 1.6, 1.6, 1.6, 1.62,
                     1.47, 1.47]),
                uncer=0.1 * SIGMA[20:38, 2],
                units='MHz/mb',
                refer='Rosenkranz, pers. comm., 2017'
            ),
            'Snr': SpectroscopicParameter(
                value=1.584e-17,
                uncer=1.584e-17 * 5 / 100,
                units='Hz*cm2/GHz2',
                refer='Pickett et al., 1998, i.e. JPL line compilation'
            ),
            'WB300': SpectroscopicParameter(
                value=0.56,
                uncer=0.05,
                units='MHz/mb',
                refer='Rosenkranz, pers. comm., 2017'
            ),
            'X11': SpectroscopicParameter(
                value=0.785,
                uncer=0.035,
                units='unitless',
                refer='Makarov et al. JQSRT 2011 -> Makarov et al. JQSRT 2008'
            ),
            'X16': SpectroscopicParameter(
                value=0.765,
                uncer=0.011,
                units='unitless',
                refer='Koshelev et al., 2016'
            ),
            'X05': SpectroscopicParameter(
                value=0.8,
                uncer=0.05,
                units='unitless',
                refer='Tretyakov et al., 2005 + Koshelev et al., 2016'
            ),
            'APU': SpectroscopicParameter(
                value=1.0,
                uncer=1.0,
                units='unitless',
                refer='Makarov et al. JQSRT 2011'
            ),
            # TODO: check better solution to referencing value within uncer
            'Y300': SpectroscopicParameter(
                value=np.array(
                    [-0.036, 0.2547, -0.3655, 0.5495, -0.5696, 0.6181, -0.4252, 0.3517, -0.1496, 0.043, 0.064, -0.1605,
                     0.2906, -0.373, 0.4169, -0.4819, 0.4963, - \
                     0.5481, 0.5512, -0.5931, 0.6212, -0.6558, 0.692, -0.7208,
                     0.7312, -0.755, 0.7555, -0.7751, 0.7914, -0.8073, 0.8307, - \
                     0.8431, 0.8676, -0.8761, 0.9046, -0.9092,
                     0.9416, -0.9423, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                uncer=1.0,
                units='1/bar == 1/1e5Pa =~ 1/atm',
                refer='Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017'
            ),
            # TODO: check better solution to referencing value within uncer
            'Y300_NL': SpectroscopicParameter(
                value=np.array(
                    [-0.036, 0.2547, -0.3655, 0.5495, -0.5696, 0.6181, -0.4252, 0.3517, -0.1496, 0.043, 0.064, -0.1605,
                     0.2906, -0.373, 0.4169, -0.4819, 0.4963, - \
                     0.5481, 0.5512, -0.5931, 0.6212, -0.6558, 0.692, -0.7208,
                     0.7312, -0.755, 0.7555, -0.7751, 0.7914, -0.8073, 0.8307, - \
                     0.8431, 0.8676, -0.8761, 0.9046, -0.9092,
                     0.9416, -0.9423, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                uncer=1.0,
                units='1/mb',
                refer='Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017'
            ),
            'O2_V': SpectroscopicParameter(
                value=np.array(
                    [0.0079, -0.0978, 0.0844, -0.1273, 0.0699, -0.0776, 0.2309, -0.2825, 0.0436, -0.0584, 0.6056, -0.6619,
                     0.6451, -0.6759, 0.6547, -0.6675, 0.6135, - \
                     0.6139, 0.2952, -0.2895, 0.2654, -0.259, 0.375, -0.368,
                     0.5085, -0.5002, 0.6206, -0.6091, 0.6526, - \
                     0.6393, 0.664, -0.6475, 0.6729, -0.6545, 0.68, -0.66,
                     0.685, -0.665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                uncer=1.0,
                units='1/bar == 1/1e5Pa =~ 1/atm',
                refer='Tretyakov et al., 2005; Uncertainty from Rosenkranz, pers. comm., 2017'
            ),
            'O2_V_NL': SpectroscopicParameter(
                value=np.array(
                    [0.0079, -0.0978, 0.0844, -0.1273, 0.0699, -0.0776, 0.2309, -0.2825, 0.0436, -0.0584, 0.6056, -0.6619,
                     0.6451, -0.6759, 0.6547, -0.6675, 0.6135, - \
                     0.6139, 0.2952, -0.2895, 0.2654, -0.259, 0.375, -0.368,
                     0.5085, -0.5002, 0.6206, -0.6091, 0.6526, - \
                     0.6393, 0.664, -0.6475, 0.6729, -0.6545, 0.68, -0.66,
                     0.685, -0.665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                uncer=1.0,
                units='1/mb',
                refer='Rosenkranz, pers. comm., 2017'
            ),
            'O2_nS': SpectroscopicParameter(
                value=2.0,
                uncer=1.0 * 0.001,
                units='unitless',
                refer='Rosenkranz, email of 2017/10/06, based on Gamache et al., 2017'
            ),
        }

        # np_zeros = np.zeros(AMU['O2gamma'].value.shape)
        # np_zeros[38:49] = AMU['O2gamma'].value[38:49] * 0.1
        # AMU['O2gamma_mmW'] = AMU['O2gamma_mmW']._replace(uncer=np_zeros)
        # np_zeros = np.zeros(AMU['O2gamma'].value.shape)
        # np_zeros[34:49] = AMU['O2gamma'].value[34:49] * 0.1
        # AMU['O2gamma_NL'] = AMU['O2gamma_NL']._replace(uncer=np_zeros)
        # np_zeros = np.zeros(AMU['O2gamma'].value.shape)
        # np_zeros[20:38] = SIGMA[20:38, 2]
        # AMU['O2gamma_WL'] = AMU['O2gamma_WL']._replace(uncer=np_zeros)
        # np_zeros = np.zeros(AMU['Y300'].value.shape)
        # np_zeros[0:34] = U[37:71]
        # AMU['Y300'] = AMU['Y300']._replace(uncer=np_zeros)
        # np_zeros = np.zeros(AMU['Y300'].value.shape)
        # np_zeros[34:49] = AMU['Y300'].value[34:49] * 0.2
        # AMU['Y300_NL'] = AMU['Y300_NL']._replace(uncer=np_zeros)
        # np_zeros = np.zeros(AMU['O2_V'].value.shape)
        # np_zeros[0:34] = U[71:105]
        # AMU['O2_V'] = AMU['O2_V']._replace(uncer=np_zeros)
        # np_zeros = np.zeros(AMU['O2_V'].value.shape)
        # np_zeros[34:49] = AMU['O2_V'].value[34:49] * 0.2
        # AMU['O2_V_NL'] = AMU['O2_V_NL']._replace(uncer=np_zeros)

        np_zeros = np.zeros(SP['O2gamma'].value.shape)
        np_zeros[38:49] = SP['O2gamma'].value[38:49] * 0.1
        SP['O2gamma_mmW'].uncer = np_zeros
        np_zeros = np.zeros(SP['O2gamma'].value.shape)
        np_zeros[34:49] = SP['O2gamma'].value[34:49] * 0.1
        SP['O2gamma_NL'].uncer = np_zeros
        np_zeros = np.zeros(SP['O2gamma'].value.shape)
        np_zeros[20:38] = SIGMA[20:38, 2]
        SP['O2gamma_WL'].uncer = np_zeros
        np_zeros = np.zeros(SP['Y300'].value.shape)
        np_zeros[0:34] = U[37:71]
        SP['Y300'].uncer = np_zeros
        np_zeros = np.zeros(SP['Y300'].value.shape)
        np_zeros[34:49] = SP['Y300'].value[34:49] * 0.2
        SP['Y300_NL'].uncer = np_zeros
        np_zeros = np.zeros(SP['O2_V'].value.shape)
        np_zeros[0:34] = U[71:105]
        SP['O2_V'].uncer = np_zeros
        np_zeros = np.zeros(SP['O2_V'].value.shape)
        np_zeros[34:49] = SP['O2_V'].value[34:49] * 0.2
        SP['O2_V_NL'].uncer = np_zeros

        if USEKOSHELEV2017:
            SP['gamma_a'] = SP['gamma_a_{}_k2017'.format(
                USEKOSHELEV2017_WHAT)]
            SP['gamma_w'] = SP['gamma_w_{}_k2017'.format(
                USEKOSHELEV2017_WHAT)]
            SP['delta_a'] = SP['delta_a_{}_k2017'.format(
                USEKOSHELEV2017_WHAT)]
            SP['delta_w'] = SP['delta_w_{}_k2017'.format(
                USEKOSHELEV2017_WHAT)]

        for i in range(0, 2):
            SP['SR'].value[i], \
                SP['SR'].uncer[i] = AbsModUncertainty.uncertainty_propagation('A/B',
                                                                              SP[
                                                                                  'delta_a'].value[i],
                                                                              SP[
                                                                                  'gamma_a'].value[i],
                                                                              SP[
                                                                                  'delta_a'].uncer[i],
                                                                              SP['gamma_a'].uncer[i])

        SP['gamma_a'].value = SP['gamma_a'].value * MB2TORR
        SP['gamma_a'].uncer = SP['gamma_a'].uncer * MB2TORR
        # AMU['gamma_a'].uncer[0:2] = np.array([0.039, 0.015])
        SP['gamma_w'].value = SP['gamma_w'].value * MB2TORR
        SP['delta_a'].value = SP['delta_a'].value * MB2TORR
        SP['delta_w'].value = SP['delta_w'].value * MB2TORR
        SP['gamma_a'].units = SP['gamma_w'].units = SP[
            'delta_a'].units = SP['delta_w'].units = 'MHz/mb'

        # units in Ros. model are [Hz*cm^2]; the conversion factor is just speed of light in cm (P. Rosenkranz, personal communication)
        # AMU['S'] = AMU['S']._replace(value=AMU['S'].value * C_CM, uncer=AMU['S'].uncer * C_CM, units='Hz*cm^2')
        SP['S'].value = SP['S'].value * C_CM
        SP['S'].uncer = SP['S'].uncer * C_CM
        SP['S'].units = 'Hz*cm^2'

        SpectroscopicParameter.SP = SP

        return SP


class AbsModUncertainty:
    """This module provides the uncertainties affecting absorption model
    coefficients found in the litterature.
    The baseline are the routines of Rosenkranz 2016 + modification to water-2-air by [Koshelev-2015]_.

    For a more exaustive example on how uncertainty work look at :ref:`uncert_example`.

    References:
        .. [1] [Koshelev-2011]_
        .. [2] [Koshelev-2015]_
        .. [3] [Koshelev-2018]_
        .. [4] [Koshelev-2017]_
        .. [5] [Turner-2009]_
        .. [6] [Tretyakov-2016]_
    """

    @staticmethod
    def uncertainty_propagation(fun: str, A: np.ndarray, B: np.ndarray, sA: np.ndarray, sB: np.ndarray,
                                a: Optional[float] = 1.0, b: Optional[float] = 1.0, sAB: Optional[float] = 0.0) -> \
            Tuple[np.ndarray, np.ndarray]:
        r"""This function propagates uncertainty given two variables A, B and their
        associated uncertainty.

        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | Function                                   |  Standard deviation                                                                                                                                                                |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | :math:`f=aA`                               |  :math:`\sigma_{f}=|a|\sigma _{A}`                                                                                                                                                 |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | :math:`\displaystyle f=aA+bB`              | :math:`\displaystyle \sigma _{f}={\sqrt {a^{2}\sigma _{A}^{2}+b^{2}\sigma _{B}^{2}+2ab\,\sigma _{AB}}}`                                                                            |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | :math:`\displaystyle f=aA-bB`              | :math:`\displaystyle \sigma _{f}={\sqrt {a^{2}\sigma _{A}^{2}+b^{2}\sigma _{B}^{2}-2ab\,\sigma _{AB}}}`                                                                            |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | :math:`\displaystyle f=AB`                 | :math:`\displaystyle \sigma _{f}\approx \left|f\right|{\sqrt {\left({\frac {\sigma _{A}}{A}}\right)^{2}+\left({\frac {\sigma _{B}}{B}}\right)^{2}+2{\frac {\sigma _{AB}}{AB}}}}`   |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | :math:`\displaystyle f={\frac {A}{B}}`     | :math:`\displaystyle \sigma _{f}\approx \left|f\right|{\sqrt {\left({\frac {\sigma _{A}}{A}}\right)^{2}+\left({\frac {\sigma _{B}}{B}}\right)^{2}-2{\frac {\sigma _{AB}}{AB}}}}`   |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
        | :math:`\displaystyle f=A^{a}`              | :math:`\displaystyle \sigma _{f}=|a|A^{a-1}\sigma _{A}`                                                                                                                            |
        +--------------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

        Args:
            fun (str): [description].
            A (numpy.ndarray): variable A
            B (numpy.ndarray): variable B
            sA (numpy.ndarray): uncertainty on A
            sB (numpy.ndarray): uncertainty on B
            a (float): multiplier for A. Defaults to 1.0.
            b (float): multiplier for B. Defaults to 1.0.
            sAB (float): sqrt of covariance between A and B. Defaults to 0.0.

        Raises:
            ValueError: Raise if `fun` argument is not defined.

        Returns:
            tuple: function computed with imput variables and uncertainty on function used

        References:
            .. [1] Wikipedia: Propagation of uncertainty https://en.wikipedia.org/wiki/Propagation_of_uncertainty
        """
        if fun == 'aA':
            F = a * A
            sF = np.abs(a) * sA
        elif fun == 'aA+bB':
            F = a * A + b * B
            sF = np.sqrt(a ** 2 * sA ** 2 + b ** 2 * sB ** 2 + 2 * a * b * sAB)
        elif fun == 'aA-bB':
            F = a * A - b * B
            sF = np.sqrt(a ** 2 * sA ** 2 + b ** 2 * sB ** 2 - 2 * a * b * sAB)
        elif fun == 'A*B':
            F = A * B
            # sF = np.sqrt( B**2*sA**2 + A**2*sB**2 + 2*A*B*sAB )
            sF = np.abs(A * B) * np.sqrt((sA / A) ** 2 +
                                         (sB / B) ** 2 + 2 * sAB / (A * B))
        elif fun == 'A/B':
            F = A / B
            sF = np.abs(A / B) * np.sqrt((sA / A) ** 2 +
                                         (sB / B) ** 2 - 2 * sAB / (A * B))
        elif fun == 'A^a':
            F = A ** a
            sF = np.abs(a) * A ** (a - 1) * sA
        else:
            raise ValueError("[utils] fun argument is mandatory")

        return F, sF

    @staticmethod
    def parameters_perturbation(what: Optional[list] = [], mode: Optional[str] = 'non', index: Optional[int] = None) -> Dict:
        """This functoin execute the perturbatoin of the spectroscopic parameters provided.

        Args:
            what (Optional[list], optional): The name of the parameter/s being perturbed. Defaults to [].
            mode (Optional[str], optional): The type of perturbation (max or min). Defaults to 'non'.
            index (Optional[int], optional): The index of the coefficient being perturbed. Defaults to None.

        Raises:
            ValueError: If argument `mode` is not set or invalid

        Returns:
            dict: The new dictionary that includes perturbed parameters.
        """
        amu_copy = deepcopy(SpectroscopicParameter.SP)
        if what == []:
            param = list(amu_copy.keys())
        else:
            param = what

        # remember here we are using matlab logig for indexing (i.e. indes start from 1 instead that 0)
        # must be changed!!!!!
        param_index = index  # - 1 if index else index

        npar = len(param)
        for i in range(npar):
            new_value = amu_copy[param[i]].value
            uncer = amu_copy[param[i]].uncer
            if mode == 'max':
                # new_value[param_index] = new_value[param_index] + uncer[param_index]
                if isinstance(new_value, float) or isinstance(uncer, float):
                    amu_copy[param[i]].value = new_value + uncer
                else:
                    amu_copy[param[i]].value[param_index] = new_value[param_index] + \
                        uncer[param_index]
            elif mode == 'min':
                # new_value[param_index] = new_value[param_index] - uncer[param_index]
                if isinstance(new_value, float) or isinstance(uncer, float):
                    amu_copy[param[i]].value = new_value - uncer
                else:
                    amu_copy[param[i]].value[param_index] = new_value[param_index] - \
                        uncer[param_index]
            elif mode == 'non':
                pass
            elif mode == 'ran':
                # here rand or randn shall be used?
                # note that rand is only positive ([0-1]) and randn gives more
                # weight to values close to zero.
                # So, one way could be to combine them, using rand and the sign of randn
                # Assuming uncorrelated uncertainty:
                # for ip1 = param_indx
                #     p1.value(ip1) = p1.value(ip1) + p1.uncer(ip1) * rand(1) * sign(randn(1));
                # end
                raise NotImplementedError(
                    "sorry, random is not implemented yet")
            else:
                raise ValueError("mode args is not set or invalid")

        return amu_copy

    # @NotImplemented
    # @staticmethod
    # def _apu_line_mixing(f: np.ndarray, t: np.ndarray, p: np.ndarray, residual_source: Optional[int] = 0):
    #     unkn = 0.0  # unknown uncertainty outside the temp and freq range
    #     if residual_source == 1:
    #         pass
    #     elif residual_source == 2:
    #         pass
    #     else:
    #         tgrid, fgrid, pgrid = np.meshgrid(TVEC, FVEC, PVEC)
    #         apu = si.interpn((tgrid, fgrid, pgrid), PR_EXT, t, f, p, 'linear')

    #     return apu

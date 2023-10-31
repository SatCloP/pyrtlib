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

from copy import deepcopy
from dataclasses import dataclass
from typing import Tuple, Optional, Dict

import numpy as np
# import scipy.interpolate as si
from pyrtlib.absorption_model import H2OAbsModel, O2AbsModel, O3AbsModel

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
        >>> from pyrtlib.uncertainty import SpectroscopicParameter
        >>> parameters = SpectroscopicParameter.water_parameters('R18')
        >>> parameters['con_Cf'].value
        5.95e-10

        New value may be added to parameters using :py:func:`~pyrtlib.uncertainty.SpectroscopicParameter` class as following

        >>> parameters['con_Xs'] = SpectroscopicParameter(2.3, 0.001, 'unitless', 'Tretyakov, JMS, 2016')
        >>> parameters['con_Xs'].value
        2.3
        >>> parameters['con_Xs'].uncer
        'unitless'
        >>> parameters['con_Xs'].refer
        'Tretyakov, JMS, 2016'

        Also, existing parameters may be modified as following

        >>> parameters['w2a'].value = 1.333
        >>> parameters['w2a'].value
        1.333

    Note:
        If new value will be added or modified it is necessary to save the new values by calling
        the :py:func:`~pyrtlib.uncertainty.SpectroscopicParameter.set_parameters` function.

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

    @staticmethod
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

        Example:
            >>> sp = SpectroscopicParameter.water_parameters("R19")
            >>> sp['con_Xs']
            SpectroscopicParameter(value=7.5, 
                                   uncer=0.0, 
                                   units='unitless', 
                                   refer='Tretyakov, JMS, 2016', 
                                   name='Self broadening temperature dependence exponents')
        """
        H2OAbsModel.model = model
        H2OAbsModel.set_ll()
        ll = H2OAbsModel.h2oll

        h2o_sp = {
            "con_Cf": SpectroscopicParameter(value=ll.cf, uncer=.0, units='1/(km*(mb^2*GHz^2))', name='Foreign induced broadening coefficient'),
            "con_Cs": SpectroscopicParameter(value=ll.cs, uncer=.0, units='1/(km*(mb^2*GHz^2))', name='Self induced broadening coefficient'),

            "con_Xf": SpectroscopicParameter(value=ll.xcf, uncer=.0, units='unitless', name='Foreign broadening temperature dependence exponents'),
            "con_Xs": SpectroscopicParameter(value=ll.xcs, uncer=.0, units='unitless', name='Self broadening temperature dependence exponents'),

            "FL": SpectroscopicParameter(value=ll.fl, uncer=np.zeros(len(ll.fl)), units='GHz', name='Central frequency line'),
            "S": SpectroscopicParameter(value=ll.s1, uncer=np.zeros(len(ll.fl)), units='Hz*cm^2', name='Line intensity (or strength)'),
            "B2": SpectroscopicParameter(value=ll.b2, uncer=.0, units='unitless', name='Temperature coefficient of intensity'),

            "gamma_a": SpectroscopicParameter(value=ll.w0*1e3, uncer=np.zeros(len(ll.fl)), units='MHz/mb', name='Air induced broadening coefficients'),
            "gamma_w": SpectroscopicParameter(value=ll.w0s*1e3, uncer=np.zeros(len(ll.fl)), units='MHz/mb', name='Water induced broadening coefficients'),
            "n_a": SpectroscopicParameter(value=ll.x, uncer=np.zeros(len(ll.fl)), units='unitless', name='Temperature-dependence exponent of air'),
            "n_w": SpectroscopicParameter(value=ll.xs, uncer=np.zeros(len(ll.fl)), units='unitless', name='Temperature-dependence exponent of water'),

            "wv_nS": SpectroscopicParameter(value=.0, uncer=.0, name='Intensity temperature-dependence exponent'),
        }
        if H2OAbsModel.model > 'R17':
            h2o_sp_ = {
                "con_Cf_factr": SpectroscopicParameter(value=1, uncer=.0, units='unitless', refer='Turner et al., TGRSS, 2009', name=''),
                "con_Cs_factr": SpectroscopicParameter(value=1, uncer=.0, units='unitless', refer='Turner et al., TGRSS, 2009', name=''),
                "n_da": SpectroscopicParameter(value=ll.xh, uncer=np.zeros(len(ll.fl)), name='T-exponent of air-shifting'),
                "n_dw": SpectroscopicParameter(value=ll.xhs, uncer=np.zeros(len(ll.fl)), name='T-exponent of water-shifting'),
                "delta_a": SpectroscopicParameter(value=ll.sh*1e3, uncer=np.zeros(len(ll.fl)), units='MHz/mb', name='Air induced shifting coefficient'),
                "delta_w": SpectroscopicParameter(value=ll.shs*1e3, uncer=np.zeros(len(ll.fl)), units='MHz/mb', name='Water induced shifting coefficient'),
            }
        else:
            h2o_sp_ = {
                "SR": SpectroscopicParameter(value=ll.sr, uncer=np.zeros(len(ll.fl)), units='unitless', name='Shift to width ratio'),
                "con_Cf_factr": SpectroscopicParameter(value=1.11, uncer=np.sqrt(0.098 ** 2 + 0.03 ** 2), units='unitless', refer='Turner et al., TGRSS, 2009', name=''),
                "con_Cs_factr": SpectroscopicParameter(value=0.79, uncer=np.sqrt(0.17 ** 2 + 0.06 ** 2), units='unitless', refer='Turner et al., TGRSS, 2009', name=''), }

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
            "O2FL": SpectroscopicParameter(value=ll.f, uncer=.0, units='GHz', name='Line frequency'),
            "O2S": SpectroscopicParameter(value=ll.s300, uncer=np.zeros(len(ll.f)), units='cm2*Hz', name='Line intensity (or strength)'),
            "O2BE": SpectroscopicParameter(value=ll.be, uncer=.0, units='unitless', name='Temperature exponent for intensity'),
            "O2gamma": SpectroscopicParameter(value=ll.w300, uncer=np.zeros(len(ll.f)), units='MHz/mb', name='Self broadening temperature dependence exponents'),
            "O2gamma_NL": SpectroscopicParameter(value=ll.w300, uncer=.0, units='MHz/mb', name='Self broadening temperature dependence exponents for neglected lines (NL)'),
            "Snr": SpectroscopicParameter(value=1.6e-17, uncer=.0, units='Hz*cm2/GHz2', name='Line intensity of non-resonant pseudo-line'),
            "WB300": SpectroscopicParameter(value=ll.wb300, uncer=.0, units='MHz/mb', name='Pressure broadening of non-resonant pseudo-line'),
            "X11": SpectroscopicParameter(value=0.785, uncer=.0, units='unitless', name='Temperature dependence of broadening coefficient'),
            "X16": SpectroscopicParameter(value=0.765, uncer=.0, units='unitless', name='Temperature dependence of broadening coefficient'),
            "X05": SpectroscopicParameter(value=ll.x, uncer=.0, units='unitless', name='Temperature dependence of broadening coefficient'),
            "O2_nS": SpectroscopicParameter(value=2.0, uncer=.0, units='unitless',  name='Intensity temperature-dependence exponent'),
            "w2a": SpectroscopicParameter(value=1.2, uncer=.0, units='unitless', name='Water-to-air broadening ratio'),
        }
        if O2AbsModel.model > 'R19':
            o2_sp_ = {
                "y0": SpectroscopicParameter(value=ll.y0, uncer=.0, name=''),
                "y1": SpectroscopicParameter(value=ll.y1, uncer=.0, name=''),
                "dnu0": SpectroscopicParameter(value=ll.dnu0, uncer=.0, name=''),
                "dnu1": SpectroscopicParameter(value=ll.dnu1, uncer=.0, name=''),
            }
        else:
            o2_sp_ = {
                "O2_V": SpectroscopicParameter(value=ll.v, uncer=np.zeros(len(ll.f)), units='1/bar == 1/1e5Pa =~ 1/atm', name='O2 line mixing temperature dependence'),
                "O2_V_NL": SpectroscopicParameter(value=ll.v, uncer=.0, units='1/mb', name='O2 line mixing temperature dependence for neglected lines (NL)'),
                "Y300": SpectroscopicParameter(value=ll.y300, uncer=np.zeros(len(ll.f)), units='1/bar == 1/1e5Pa =~ 1/atm', name='Line mixing coefficients (single lines)'),
                "Y300_NL": SpectroscopicParameter(value=ll.y300, uncer=.0, units='1/mb', name='Line mixing coefficients for neglected lines (NL)'),
            }
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
            "O3_FL": SpectroscopicParameter(value=ll.fl, uncer=np.zeros(len(ll.fl)), units='GHz', name='Line frequency'),
            "O3_S1": SpectroscopicParameter(value=ll.s1, uncer=np.zeros(len(ll.fl)), units='Hz*cm^2', name='Line intensity (or strength)'),
            "O3_B": SpectroscopicParameter(value=ll.b, uncer=np.zeros(len(ll.fl)), units='unitless', name='Temperature Coefficient of intensity'),
            "O3_W": SpectroscopicParameter(value=ll.w, uncer=np.zeros(len(ll.fl)), units='GHz/mb', name='Air-pressure broadening'),
            "O3_X": SpectroscopicParameter(value=ll.x, uncer=np.zeros(len(ll.fl)), units='unitless', name='Temperature exponent of air broadening'),
            "O3_SR": SpectroscopicParameter(value=ll.sr, uncer=np.zeros(len(ll.fl)), units='unitless', name='Shift-to-width ratio'),
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
        try:
            amu_copy = deepcopy(SpectroscopicParameter.SP)
        except AttributeError as e:
            print('Spectroscopic parameters are not defined')
            raise e
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

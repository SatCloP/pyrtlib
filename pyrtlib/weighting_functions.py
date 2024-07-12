# -*- coding: utf-8 -*-
"""
This module is  used to define the weighting functions
"""

__author__ = ''
__date__ = 'May 2024'
__copyright__ = '(C) 2024, CNR-IMAA'

from typing import Dict, Tuple, Optional, List
import warnings
# warnings.filterwarnings("ignore", category=UserWarning)
import numpy as np
import matplotlib.pyplot as plt
from . absorption_model import AbsModel
from .tb_spectrum import TbCloudRTE
from .absorption_model import O2AbsModel


class WeightingFunctions(object):
    """
    This class is used to compute the weighting functions

    .. note::
        The weighting functions are computed always using last absorption model available.
    """

    def __init__(self, z: np.ndarray, p: np.ndarray, t: np.ndarray, rh: np.ndarray,
                 value_prof_interpolation: Optional[float] = 1) -> None:
        """
        Constructor for the WeightingFunctions class.

        Args:
            z (np.ndarray): Height profile.
            p (np.ndarray): Pressure profile.
            t (np.ndarray): Temperature profile.
            rh (np.ndarray): Relative humidity profile.
            value_prof_interpolation (Optional[float], optional): The value for the interpolation of the profiles. Defaults to 1.

        Raises:
            ValueError: If the input values are not valid arrays.

        Returns:
            None

        .. note::
            The input arrays must have the same size. Units for the arrays are:
                - z: km
                - p: hPa
                - t: K
                - rh: fraction
                - value_prof_interpolation: km
        """

        self.zinterp = np.arange(
            0, np.max(z)+value_prof_interpolation, value_prof_interpolation)
        self.pinterp = np.interp(self.zinterp, z, p)
        self.tinterp = np.interp(self.zinterp, z, t)
        self.rhinterp = np.interp(self.zinterp, z, rh)

        self._frequencies = None
        self._satellite = True
        self.emissivity = 1.
        self.model = AbsModel.implemented_models()['WaterVapour'][-1]
        self.angle = 90.
        self.bandpass = None
        self.legend_labels = None

    @property
    def satellite(self) -> bool:
        """If :code:`True` computes an upward-propagating brightness-temperature spectrum
        otherwise a downward-propagating brightness-temperature 
        spectrum at the bottom of the atmosphere will be performed.
        """
        return self._satellite

    @satellite.setter
    def satellite(self, sat: bool) -> None:
        """
        Sets the satellite flag.

        Parameters:
        - sat (bool): A boolean value indicating whether the satellite flag should be set to True or False.

        Raises:
        - ValueError: If the input value for `sat` is not a boolean.

        Returns:
        - None
        """
        if isinstance(sat, bool):
            self._satellite = sat
        else:
            raise ValueError("Please enter True or False")

    @property
    def frequencies(self) -> np.ndarray:
        """Returns the frequencies array.

        Returns:
            np.ndarray: The frequencies array.
        """
        return self._frequencies

    @frequencies.setter
    def frequencies(self, frequencies: np.ndarray) -> None:
        """
        Set the frequencies for the weighting function.

        Parameters:
        - frequencies: np.ndarray
            The frequencies to be set.

        Raises:
        - ValueError: If the input is not a valid array for frequencies.

        Returns:
        - None
        """
        if isinstance(frequencies, np.ndarray):
            self._frequencies = frequencies
        else:
            raise ValueError(
                "Please enter a valid array for frequencies")

    def _perform_calculation(self) -> Dict[str, np.ndarray]:
        """
        Perform the calculation to get opacity and absorption.
        """

        rte = TbCloudRTE(self.zinterp, self.pinterp, self.tinterp,
                         self.rhinterp, self.frequencies, angles=np.array([self.angle]), ray_tracing=True)
        rte.satellite = self.satellite
        rte.init_absmdl(self.model)
        rte.emissivity = self.emissivity
        O2AbsModel.model = AbsModel.implemented_models()['Oxygen'][-1]
        O2AbsModel.set_ll()
        _, other = rte.execute(False)

        return other

    def _mean_channel(self, a: np.ndarray, bandpass: np.ndarray) -> np.ndarray:
        """
        Compute the mean of the input array based on bandpass.

        Args:
            a (np.ndarray): The input array.
            bandpass (np.ndarray): The bandpass array.

        Returns:
            np.ndarray: The mean of the input array based on bandpass.
        """

        channel_mean = np.zeros((len(bandpass), a.shape[1]))
        cnt = 0
        for i, _ in enumerate(bandpass):
            channel_mean[i, :] = np.mean(a[cnt:bandpass[i]+cnt, :], axis=0)
            cnt += bandpass[i]

        return channel_mean

    def generate_wf(self) -> np.ndarray:
        r"""
        Generate the weighting function.

        .. math:: \frac{\partial \tau}{\partial s} = -{\tau (\nu, s)}\times{\alpha (\nu, s)}

        This method calculates the weighting function based on the provided data.
        It performs calculations and returns the resulting weighting function.

        Returns:
            np.ndarray: The calculated weighting function.

        """
        other = self._perform_calculation()
        if self.satellite:
            z0 = self.zinterp[-1]
        else:
            z0 = self.zinterp[0]
        self.l0 = np.argmin(np.abs(self.zinterp-z0))

        tautotal = (other['taulaywet'] + other['taulaydry'] +
                    other['taulayliq'] + other['taulayice'])
        abstotal = other['awet'] + other['adry']

        if self.satellite:
            ao = np.fliplr(tautotal[:, 0, :self.l0+1])
            abss = abstotal[:, 0, :self.l0+1]
        else:
            ao = tautotal[:, 0, self.l0:]
            abss = abstotal[:, 0, self.l0:]

        if self.satellite:
            cumao = np.fliplr(np.cumsum(ao[::-1], axis=1)[::-1])
        else:
            cumao = np.cumsum(ao, axis=1)
        trans = np.exp(-cumao)
        wgt = np.array([trans[:, j] * abss[:, j]
                       for j in range(trans.shape[1])]).T

        return wgt

    def plot_wf(self, wgt: np.ndarray, title: Optional[str] = '', ylim: Optional[Tuple[int, int]] = None,
                xlim: Optional[Tuple[int, int]] = None,
                legend: Optional[bool] = True,
                normalized: Optional[bool] = False,
                **kwargs) -> None:
        """
        Plot the weighting functions

        Args:
            wgt (ndarray): The weighting functions to be plotted.
            title (Optional[str], optional): The title of the plot. Defaults to ''.
            ylim (Optional[Tuple[int, int]], optional): The y-axis limits of the plot. Defaults to None.
            xlim (Optional[Tuple[int, int]], optional): The x-axis limits of the plot. Defaults to None.
            legend (Optional[bool], optional): Whether to show the legend. Defaults to True.
            normalized (Optional[bool], optional): Whether to normalize to the max the weighting functions. Defaults to False.
            **kwargs: Additional keyword arguments (figsize, dpi)

        """
        # plt.rcdefaults()
        plt.rcParams.update({'font.size': 15})
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['font.stretch'] = 'condensed'
        # plt.rcParams["font.weight"] = "bold"
        # plt.rcParams["axes.labelweight"] = "bold"
        # plt.rc('font', family=['cmr10'])
        # plt.rc('mathtext', fontset='cm')
        # plt.rc('axes.formatter', use_mathtext=False)

        if self.bandpass is not None:
            wgt = self._mean_channel(wgt, self.bandpass)
            self.frequencies = np.array(
                [np.mean(self.frequencies[i:i+j]) for i, j in enumerate(self.bandpass)])

        # if yaxis == 'p':
        #     y = self.pinterp[:-1]
        #     y_label = 'Pressure [$hPa$]'
        #     units = 'hPa'
        # else:
        y = self.zinterp[:-1]
        y_label = 'Height [$km$]'
        units = 'km'

        normalize = np.max(
            wgt, axis=1) if normalized else np.ones(wgt.shape[1])

        plt.figure(figsize=kwargs.get('figsize', (5, 10)),
                   dpi=kwargs.get('dpi', 150))
        for i, frequenciy in enumerate(self.frequencies):
            if self.satellite:
                line_s = 'dashed' if y[np.argmax(wgt[i, 1:])] == 0. else None
                plt.plot(wgt[i, 1:]/normalize[i], y, linestyle=line_s,
                         label=f'{frequenciy:.2f} GHz {y[np.argmax(wgt[i, 1:])]:.0f} {units}')
            else:
                plt.plot(wgt[i, :], self.zinterp,
                         label=f'{frequenciy:.2f} GHz')

        # if yaxis == 'p':
        #     plt.gca().invert_yaxis()
        #     plt.yscale('log')
        plt.xlabel('Weighting Function [$km^{-1}$]')
        plt.ylabel(y_label)
        if legend:
            if self.legend_labels:
                plt.legend(self.legend_labels, fontsize=10)
            else:
                plt.legend(fontsize=10)
        if ylim:
            plt.ylim(ylim[0], ylim[1])
        if xlim:
            plt.xlim(xlim[0], xlim[1])

        if title:
            plt.title(title)
        plt.grid(color='#d6d6d6', linestyle='dashed')
        plt.tick_params(axis='both', which='major', direction='in')
        ax_twinx = plt.twinx()
        ax_twinx.set_ylim((self.pinterp[np.where(self.zinterp == np.min(self.zinterp))],
                           self.pinterp[np.where(self.zinterp == np.max(self.zinterp))]))
        if ylim:
            ax_twinx.set_ylim((self.pinterp[np.where(self.zinterp == ylim[0])],
                               self.pinterp[np.where(self.zinterp == ylim[1])]))
        ax_twinx.set_yscale('log')
        ax_twinx.set_ylabel('Pressure [$hPa$]')
        ax_twinx.tick_params(axis='both', which='both', direction='in')
        plt.tight_layout()
        plt.show()

    def plot_wf_grouped(self, wgt: np.ndarray, title: Optional[str] = '',
                        ylim: Optional[Tuple[int, int]] = None,
                        grouped_frequencies: List[int] = None,
                        grouped_labels: List[str] = None,
                        **kwargs) -> None:
        """
        Plot the weighting functions

        Args:
            wgt (ndarray): The weighting functions to be plotted.
            title (Optional[str], optional): The title of the plot. Defaults to ''.
            ylim (Optional[Tuple[int, int]], optional): The y-axis limits of the plot. Defaults to None.
            grouped_frequencies (List[int]): The number of frequencies to be grouped.
            grouped_labels (List[str]): The labels for the grouped frequencies.
            dpi (Optional[int], optional): The dpi of the plot. Defaults to 150.
            **kwargs: Additional keyword arguments (figsize, dpi)

        Raises:
            ValueError: If the plot is not available for the current satellite mode.
        """
        if not self.satellite:
            raise ValueError("This plot is only available for satellite mode")

        # plt.rcdefaults()
        plt.rcParams.update({'font.size': 15})
        # plt.rc('font', family=['Sans Serif'])
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['font.stretch'] = 'condensed'
        # plt.rcParams["font.weight"] = "bold"
        # plt.rcParams["axes.labelweight"] = "bold"
        # plt.rc('mathtext', fontset='cm')
        # plt.rc('axes.formatter', use_mathtext=False)

        y = self.zinterp[:-1]
        y_label = 'Height [$km$]'

        plt.figure(figsize=kwargs.get('figsize', (10, 4)),
                   dpi=kwargs.get('dpi', 150))
        cnt = 0
        max_increment = np.array([])
        for i in grouped_frequencies:
            a = np.tile(cnt, i)
            max_increment = np.append(max_increment, a)
            cnt += 0.40

        cnt = 0
        g_labels = np.array([], dtype='object')
        for i, grouped_frequency in enumerate(grouped_frequencies):
            a = np.repeat(grouped_labels[i], grouped_frequency)
            g_labels = np.append(g_labels, a)
            cnt += 1

        for i in range(len(self.frequencies)):
            line_s = 'dashed' if y[np.argmax(wgt[i, 1:])] == 0. else None
            plt.plot(wgt[i, 1:]+max_increment[i], y,
                     linestyle=line_s, linewidth=1.5)

        for i in range(len(grouped_frequencies)):
            index_wgt = np.cumsum(np.array(grouped_frequencies)) - 1
            height = ylim[1]-1 if ylim else np.max(y)-1
            perc_height = height/15
            incr_wgt = np.min(wgt[index_wgt[i], 1:]) + \
                max_increment[index_wgt[i]]
            plt.text(.02 + incr_wgt, height-perc_height,
                     f'{grouped_labels[i]}', ha='left', va='top', fontsize=12)

        plt.xlabel('Weighting Function [$km^{-1}$]')
        plt.ylabel(y_label)
        if ylim:
            plt.ylim(ylim[0], ylim[1])
        if title:
            plt.title(title)
        plt.grid(color='#d6d6d6')
        plt.tick_params(axis='both', which='major', direction='in')
        ax_twinx = plt.twinx()
        ax_twinx.set_ylim((self.pinterp[np.where(self.zinterp == np.min(self.zinterp))],
                           self.pinterp[np.where(self.zinterp == np.max(self.zinterp))]))
        if ylim:
            ax_twinx.set_ylim((self.pinterp[np.where(self.zinterp == ylim[0])],
                               self.pinterp[np.where(self.zinterp == ylim[1])]))
        ax_twinx.set_yscale('log')
        ax_twinx.set_ylabel('Pressure [$hPa$]')
        ax_twinx.tick_params(axis='both', which='both', direction='in')
        plt.tight_layout()
        plt.show()

        warnings.warn(
            "Weighting functions have been shifted to the right so that they all can be displayed in a single graph")

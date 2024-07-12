from unittest import TestCase

import numpy as np
from pyrtlib.weighting_functions import WeightingFunctions
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.utils import ppmv2gkg, mr2rh
import numpy as np
import matplotlib
matplotlib.use('Template')

class Test(TestCase):
    def wf_computation(self, frq: np.ndarray):
        z, p, d, t, md = atmp.gl_atm(atmp.TROPICAL)
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh = mr2rh(p, t, gkg)[0] / 100

        wf = WeightingFunctions(z, p, t, rh)
        wf.frequencies = frq

        return wf.generate_wf(), z, wf

    def test_f_89(self):
        wgt, z, _ = self.wf_computation(np.array([89.0]))
        np.testing.assert_almost_equal(z[np.argmax(wgt)], 0., decimal=15)

    def test_f_57(self):
        wgt, z, _ = self.wf_computation(np.array([57.290344]))
        np.testing.assert_almost_equal(z[np.argmax(wgt)], 18., decimal=15)

    def test_f_57_6(self):
        wgt, z, _ = self.wf_computation(np.array([57.660544]))
        np.testing.assert_almost_equal(z[np.argmax(wgt)], 27.5, decimal=15)

    def test_plot_wf(self):
        wgt, z, wf = self.wf_computation(np.array([57.660544]))
        wf.plot_wf(wgt, 'Test', legend=True, ylim=(0, 60), xlim=(0, .1), figsize=(8, 6), dpi=100)

    def test_plot_wf_grouped(self):
        wgt, z, wf = self.wf_computation(np.array([57.660544]))
        wf.plot_wf_grouped(wgt, 'Test', grouped_frequencies=[
                           1], grouped_labels=['Test'], dpi=100)

    def test_plot_wf_multiple(self):
        freqs = np.array([57.290344, 57.660544, 89.0])
        wgt, z, wf = self.wf_computation(freqs)
        wf.satellite = False
        wf.plot_wf(wgt, 'Multiple Frequencies',
                   legend=True, figsize=(8, 6), dpi=100)
        
    def test_plot_wf_bandpass(self):
        cf53 = 53.596
        cf57 = 57.290344
        freqs = np.array([52.8, cf53-0.115, cf53+0.115, 54.4, 54.94, 55.5, cf57, cf57-0.217, cf57+0.217, cf57-0.3222-0.048, cf57-0.3222+0.048, cf57+0.3222-0.048, cf57+0.3222+0.048])
        wgt, z, wf = self.wf_computation(freqs)
        wf.bandpass = np.array([1, 2, 1, 1, 1, 1, 2, 4])
        wf.plot_wf(wgt, 'Bandpass',
                   legend=True, figsize=(8, 6), dpi=100)

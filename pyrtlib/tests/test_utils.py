from unittest import TestCase

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.utils import (ppmv2gkg, mr2rh, gas_mass, height_to_pressure, pressure_to_height,
                           to_kelvin, to_celsius, get_frequencies, eswat_goffgratch, satvap, satmix)

z, p, d, tk, md = atmp.gl_atm(atmp.TROPICAL)


class Test(TestCase):
    def test_gas_mass(self):
        # H2O
        mass_proton = 1.6726485e-27
        mass_molecule = np.dot(mass_proton, np.dot(2, 1) + 16)
        h20 = gas_mass(atmp.H2O)
        assert_allclose(h20, mass_molecule, atol=0)
        # CO2
        mass_molecule = np.dot(mass_proton, 12 + np.dot(2, 16))
        co2 = gas_mass(atmp.CO2)
        assert_allclose(co2, mass_molecule, atol=0)

    def test_mr2rh(self):
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rh, rh_wmo = mr2rh(p, tk, gkg)
        rh_ex = np.array([7.37904947e-01, 7.15135305e-01, 7.35076661e-01, 4.79160577e-01,
                          3.48167127e-01, 3.76644483e-01, 3.48042202e-01, 3.20269720e-01,
                          2.95282900e-01, 2.54125093e-01, 1.95488576e-01, 1.31680781e-01,
                          9.25731286e-02, 5.87840339e-02, 7.40605521e-02, 9.93775545e-02,
                          1.69111933e-01, 1.94814864e-01, 8.37051182e-02, 3.76010914e-02,
                          1.81872123e-02, 9.21185036e-03, 5.02984723e-03, 3.31147808e-03,
                          2.40076896e-03, 1.61089247e-03, 6.48478175e-04, 2.82194238e-04,
                          1.22168272e-04, 5.45535815e-05, 2.50874197e-05, 1.17670685e-05,
                          5.71681338e-06, 2.78356544e-06, 1.45692733e-06, 1.04299888e-06,
                          9.37799033e-07, 1.14946062e-06, 2.58322551e-06, 6.76562179e-06,
                          2.12214737e-05, 9.12535389e-05, 1.01732056e-04, 2.65502994e-05,
                          1.60922893e-06, 1.62038182e-07, 2.70297183e-09, 4.12976939e-11,
                          2.49011743e-13, 3.49656489e-15])
        assert_allclose(rh / 100, rh_ex, atol=0.001)
        rh_wmo_ex = np.array([7.31108823e-01, 7.09583292e-01, 7.31012737e-01, 4.74681358e-01,
                              3.45272337e-01, 3.74558736e-01, 3.46672438e-01, 3.19393547e-01,
                              2.94744708e-01, 2.53819434e-01, 1.95334753e-01, 1.31617341e-01,
                              9.25467678e-02, 5.87747158e-02, 7.40547927e-02, 9.93739520e-02,
                              1.69109441e-01, 1.94812529e-01, 8.37025984e-02, 3.75985892e-02,
                              1.81846596e-02, 9.20922477e-03, 5.02706131e-03, 3.30858768e-03,
                              2.39757664e-03, 1.60764771e-03, 6.44880510e-04, 2.78195367e-04,
                              1.17868797e-04, 4.99538324e-05, 2.01875426e-05, 6.56712965e-06,
                              2.16844820e-07, -2.91641869e-06, -4.44306407e-06, -4.95699486e-06,
                              -5.06219534e-06, -4.85053248e-06, -2.81676054e-06, 2.26565223e-06,
                              1.79215437e-05, 8.91537306e-05, 1.00432188e-04, 2.57003220e-05,
                              1.06922980e-06, -2.37961753e-07, -3.37297027e-07, -2.79958702e-07,
                              -2.39999751e-07, -1.99999997e-07])
        assert_allclose(rh_wmo / 100, rh_wmo_ex, atol=0)

    def test_ppmv2gkg(self):
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        gkg_ex = np.array([1.61276973e+01, 1.21222067e+01, 9.54102877e+00, 5.34894703e+00,
                           2.76217137e+00, 2.08111358e+00, 1.30676020e+00, 8.01720084e-01,
                           4.74998936e-01, 2.54883546e-01, 1.18920776e-01, 4.54411709e-02,
                           1.80682455e-02, 6.15750879e-03, 3.86865704e-03, 2.48788234e-03,
                           1.86591175e-03, 1.80371470e-03, 1.71041911e-03, 1.61712352e-03,
                           1.61712352e-03, 1.64822205e-03, 1.74151764e-03, 1.80371470e-03,
                           1.99030587e-03, 2.02140440e-03, 2.23909411e-03, 2.48788234e-03,
                           2.67447352e-03, 2.86106469e-03, 3.04765587e-03, 3.23424704e-03,
                           3.42083822e-03, 3.54523233e-03, 3.66962645e-03, 3.73182351e-03,
                           3.73182351e-03, 3.73182351e-03, 3.35864116e-03, 2.79886763e-03,
                           2.05250293e-03, 1.30613823e-03, 8.08561760e-04, 5.28674997e-04,
                           3.35864116e-04, 2.48788234e-04, 2.11469999e-04, 1.74151764e-04,
                           1.49272940e-04, 1.24394117e-04])
        assert_allclose(gkg, gkg, atol=0)

    def test_pressure_to_height(self):
        h = pressure_to_height(850)
        assert_almost_equal(h, 1506.02, decimal=2)

    def test_height_to_pressure(self):
        p = height_to_pressure(5000)
        assert_almost_equal(p, 577.63, decimal=2)

    def test_to_kelvin(self):
        t = to_kelvin(25)
        assert_almost_equal(t, 298.25, decimal=2)

    def test_to_celsius(self):
        t = to_celsius(298.25)
        assert_almost_equal(t, 25, decimal=2)

    def test_get_frequencies(self):
        frequencies = np.array([18.7,  23.8,  31.4,  50.3,  52.61,  53.24,
                                53.75,  89., 115.5503, 116.6503, 117.3503, 117.5503,
                                119.9503, 120.1503, 120.8503, 121.9503, 164.775, 166.225,
                                174.91, 177.21, 178.41, 179.91, 181.31, 185.31,
                                186.71, 188.21, 189.41, 191.71])
        
        assert_almost_equal(frequencies, get_frequencies('MWI'), decimal=4)

    def test_eswat_goffgratch(self):
        eswat = eswat_goffgratch(273.25)
        assert_almost_equal(eswat, 6.147856307200233, decimal=10)
        eswat = eswat_goffgratch(np.array([273.25, 280.5]))
        assert_almost_equal(eswat, np.array([ 6.14785631, 10.24931596]), decimal=5)

    def test_satvap(self):
        sat_vap = satvap(273.25)
        assert_almost_equal(sat_vap, 6.147856307200233, decimal=10)
        sat_vap = satvap(np.array([273.25, 280.5]))
        assert_almost_equal(sat_vap, np.array([ 6.14785631, 10.24931596]), decimal=5)

    def test_satmix(self):
        sat_mix = satmix(1000, 273.25)
        assert_almost_equal(sat_mix, 3.8474392877771995, decimal=10)
        sat_mix = satmix(np.array([1000, 950, 850]), np.array([273.25, 260.34, 258.36]))
        assert_almost_equal(sat_mix, np.array([3.84743929, 1.4993625 , 1.42541609]), decimal=5)

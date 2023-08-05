from unittest import TestCase

import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal, assert_equal
from pyrtlib.absorption_model import H2OAbsModel
from pyrtlib.climatology import AtmosphericProfiles as atmp
from pyrtlib.utils import (ppmv2gkg, mr2rh, gas_mass, height_to_pressure, pressure_to_height, constants,
                           to_kelvin, to_celsius, get_frequencies_sat, eswat_goffgratch, satvap, satmix,
                           import_lineshape, atmospheric_tickness, mr2rho, mr2e, esice_goffgratch)

z, p, d, t, md = atmp.gl_atm(atmp.TROPICAL)


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
        rh, rh_wmo = mr2rh(p, t, gkg)
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

    def test_mr2rho(self):
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        rho = mr2rho(gkg, t, p)
        rho_ex = np.array([1.85104494e+01, 1.27497840e+01, 9.15966006e+00, 4.65624332e+00,
                           2.18922176e+00, 1.49434354e+00, 8.47898228e-01, 4.68871213e-01,
                           2.49707548e-01, 1.19873202e-01, 4.99841165e-02, 1.69917359e-02,
                           5.99585055e-03, 1.79908641e-03, 9.99729444e-04, 5.61630430e-04,
                           3.66257206e-04, 3.02243417e-04, 2.36483910e-04, 1.85098191e-04,
                           1.53988988e-04, 1.30806974e-04, 1.15627166e-04, 1.01347846e-04,
                           9.48940031e-05, 8.17423145e-05, 6.05810598e-05, 4.55174946e-05,
                           3.33953866e-05, 2.45998013e-05, 1.82006388e-05, 1.35293422e-05,
                           1.01070051e-05, 7.41585104e-06, 5.50043485e-06, 4.10894800e-06,
                           2.25064600e-06, 1.22761977e-06, 5.99893638e-07, 2.58346312e-07,
                           9.21242928e-08, 2.70843598e-08, 6.99820044e-09, 1.78971225e-09,
                           4.36783378e-10, 1.31345876e-10, 4.51747476e-11, 1.62470581e-11,
                           6.24649461e-12, 2.56589160e-12])
        assert_allclose(rho, rho_ex, atol=0.001)

    def test_mr2e(self):
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        e = mr2e(p, gkg)
        e_ex = np.array([2.56031990e+01, 1.72821313e+01, 1.21621329e+01, 6.09656950e+00,
                         2.79872387e+00, 1.86417647e+00, 1.03152477e+00, 5.56131147e-01,
                         2.88458304e-01, 1.34768972e-01, 5.46727466e-02, 1.80445017e-02,
                         6.18747025e-03, 1.80178216e-03, 9.70313965e-04, 5.27997888e-04,
                         3.32999001e-04, 2.71729212e-04, 2.16974403e-04, 1.73159550e-04,
                         1.46899618e-04, 1.27199663e-04, 1.14519679e-04, 1.01499706e-04,
                         9.59996928e-05, 8.35247285e-05, 6.34677715e-05, 4.87998048e-05,
                         3.66358425e-05, 2.75998730e-05, 2.08738977e-05, 1.58599175e-05,
                         1.20999335e-05, 9.06294834e-06, 6.84395962e-06, 5.12396926e-06,
                         2.73598358e-06, 1.43399140e-06, 6.53396472e-07, 2.60998826e-07,
                         8.57997169e-08, 2.30999515e-08, 5.71999256e-09, 1.46199876e-09,
                         3.71519799e-10, 1.15599954e-10, 4.41999850e-11, 1.81159949e-11,
                         8.63999793e-12, 4.49999910e-12])
        assert_allclose(e, e_ex, atol=0.001)

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

    def test_get_frequencies_sat(self):
        frequencies = np.array([18.7,  23.8,  31.4,  50.3,  52.61,  53.24,
                                53.75,  89., 115.5503, 116.6503, 117.3503, 117.5503,
                                119.9503, 120.1503, 120.8503, 121.9503, 164.775, 166.225,
                                174.91, 177.21, 178.41, 179.91, 181.31, 185.31,
                                186.71, 188.21, 189.41, 191.71])

        assert_almost_equal(frequencies, get_frequencies_sat('MWI'), decimal=4)

    def test_eswat_goffgratch(self):
        eswat = eswat_goffgratch(273.25)
        assert_almost_equal(eswat, 6.147856307200233, decimal=10)
        eswat = eswat_goffgratch(np.array([273.25, 280.5]))
        assert_almost_equal(eswat, np.array(
            [6.14785631, 10.24931596]), decimal=5)

    def test_satvap(self):
        sat_vap = satvap(273.25)
        assert_almost_equal(sat_vap, 6.147856307200233, decimal=10)
        sat_vap = satvap(np.array([273.25, 280.5]))
        assert_almost_equal(sat_vap, np.array(
            [6.14785631, 10.24931596]), decimal=5)

    def test_satmix(self):
        sat_mix = satmix(1000, 273.25)
        assert_almost_equal(sat_mix, 3.8474392877771995, decimal=10)
        sat_mix = satmix(np.array([1000, 950, 850]),
                         np.array([273.25, 260.34, 258.36]))
        assert_almost_equal(sat_mix, np.array(
            [3.84743929, 1.4993625, 1.42541609]), decimal=5)

    def test_constants(self):
        cs = constants('R')[0]
        assert_almost_equal(cs, 8.31446261815324, decimal=5)
        cs_list = constants()
        cs_list_ex = ['avogadro',
                      'boltzmann',
                      'EarthRadius',
                      'gravity',
                      'light',
                      'Np2dB',
                      'planck',
                      'Rdry',
                      'Rwatvap',
                      'Tcosmicbkg',
                      'R']
        assert_equal(cs_list, cs_list_ex)

    def test_import_linelist(self):
        model = 'R21SD'
        h2oll = import_lineshape('h2oll_{}'.format(model))
        aself = np.array([0., 12.6,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                          0.,  0.,  0.,  0.,  0.])
        assert_almost_equal(h2oll.aself, aself, decimal=5)
        H2OAbsModel.model = 'R21SD'
        H2OAbsModel.set_ll()
        assert_almost_equal(H2OAbsModel.h2oll.aself, aself, decimal=5)

    def test_atmospheric_tickness(self):
        z = atmospheric_tickness(np.array([1000, 500]), np.array([280, 244]))
        assert_almost_equal(np.array([0., 5.31569228]), z, decimal=5)
        gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
        z = atmospheric_tickness(p, t, gkg/1000)
        h = np.array([0.,   0.99705536,   1.9904238,   2.98635587,
                      3.9884136,   4.98568001,   5.98430371,   6.97583587,
                      7.96762711,   8.97141576,   9.95671685,  10.95897373,
                      11.94235997,  12.95659692,  13.9206096,  14.93279803,
                      15.94893135,  16.92047757,  17.91080171,  18.90666321,
                      19.89208472,  20.88806132,  21.88441018,  22.86842643,
                      23.85251614,  24.85011595,  27.3234496,  29.79827113,
                      32.26783235,  34.73529522,  37.19942755,  39.65669422,
                      42.1113353,  44.60256935,  47.06868679,  49.48806998,
                      54.38801236,  59.27146052,  64.14382998,  69.03950383,
                      73.97963167,  78.84667776,  83.69984482,  88.56755124,
                      93.4126703,  98.17297656, 102.88137249, 107.51360655,
                      112.1578853, 116.83331575])
        assert_almost_equal(h, z, decimal=5)

    def test_esice_goffgratch(self):
        eice = esice_goffgratch(t)
        eice_ex = np.array([4.46886926e+01, 2.94273085e+01, 1.90367102e+01, 1.40918238e+01,
                            8.34828341e+00, 4.81745835e+00, 2.70260886e+00, 1.48458452e+00,
                            7.82312951e-01, 3.97932076e-01, 1.96960374e-01, 9.04389659e-02,
                            4.15819757e-02, 1.80168867e-02, 7.30775615e-03, 2.83619102e-03,
                            1.01755016e-03, 7.15761625e-04, 1.34919490e-03, 2.44422502e-03,
                            4.39341246e-03, 7.72458595e-03, 1.31247892e-02, 1.80168867e-02,
                            2.39426078e-02, 3.16390238e-02, 6.27796299e-02, 1.16494842e-01,
                            2.12606526e-01, 3.77790389e-01, 6.54739552e-01, 1.11904734e+00,
                            1.85249073e+00, 3.00390376e+00, 4.54128870e+00, 4.77709067e+00,
                            2.65516657e+00, 1.02671926e+00, 1.76446210e-01, 2.30397076e-02,
                            2.13532575e-03, 1.30308709e-04, 3.08419643e-05, 3.02455186e-05,
                            1.19097392e-04, 3.63686345e-04, 9.23737220e-03, 3.22861540e-01,
                            4.46886926e+01, 3.25249157e+03])
        assert_almost_equal(eice_ex, eice, decimal=5)

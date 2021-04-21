from unittest import TestCase

from numpy.testing import assert_allclose
from pyrtlib.absmod_uncertainty import absmod_uncertainties_perturb


class Test(TestCase):
    def test_absmod_uncertainties_perturb_min(self):
        amu = absmod_uncertainties_perturb(['gamma_a'], 'min', 2)
        expected0 = 0.016501364000000
        assert_allclose(expected0, amu['gamma_a'].uncer[0])
        expected1 = 0.015001240000000
        assert_allclose(expected1, amu['gamma_a'].uncer[1])

        expected0 = 2.698723076000000
        assert_allclose(expected0, amu['gamma_a'].value[0])
        expected1 = 2.929742172000000
        assert_allclose(expected1, amu['gamma_a'].value[1])

    def test_absmod_uncertainties_perturb_max(self):
        amu = absmod_uncertainties_perturb(['gamma_a'], 'max', 1)
        expected0 = 2.715224440000000
        assert_allclose(expected0, amu['gamma_a'].value[0])

        amu = absmod_uncertainties_perturb(['gamma_a'], 'max', 2)
        expected1 = 2.959744652000000
        assert_allclose(expected1, amu['gamma_a'].value[1])

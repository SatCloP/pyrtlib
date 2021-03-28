from unittest import TestCase
import numpy as np
from numpy.testing import assert_allclose

from pyrtlib.utils import gas_mass

class Test(TestCase):
    def test_gas_mass(self):
        # H2O
        mass_proton = 1.6726485e-27
        mass_molecule = np.dot(mass_proton, np.dot(2, 1) + 16)
        h20 = gas_mass(1)
        assert_allclose(h20, mass_molecule, atol=0.001)
        # CO2
        mass_molecule = np.dot(mass_proton, 12 + np.dot(2, 16))
        co2 = gas_mass(2)
        assert_allclose(co2, mass_molecule, atol=0.001)


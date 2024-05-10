import os
# from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose

from unittest import TestCase
from pyrtlib.absorption_model import (H2OAbsModel, O2AbsModel,
                                      N2AbsModel, LiqAbsModel)

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class Test(TestCase):
    def test_setter_model(self):
        h2o = H2OAbsModel()
        h2o.model = 'R22SD'
        assert 'R22' != h2o.model

        o2 = O2AbsModel()
        o2.model = 'R22'
        assert 'R20' != o2.model

    def test_lineshape(self):
        H2OAbsModel.model = 'R22SD'
        H2OAbsModel.set_ll()
        assert H2OAbsModel.h2oll.reftline == 296.

        H2OAbsModel.model = 'R98'
        H2OAbsModel.set_ll()
        assert H2OAbsModel.h2oll.reftline != 296.

        if hasattr(H2OAbsModel, 'model'):
            assert H2OAbsModel.model != 'PIPPO'

    def test_absn2(self):
        N2AbsModel.model = 'R22SD'
        absn2 = N2AbsModel.n2_absorption(267., 800., 55.0034)
        assert absn2 == 0.000278315910216229

        N2AbsModel.model = 'R98'
        absn2 = N2AbsModel.n2_absorption(267., 800., 55.0034)
        assert absn2 != 0.000278315910216229

    def test_absliq(self):
        LiqAbsModel.model = 'R22SD'
        absliq = LiqAbsModel.liquid_water_absorption(0.05, 183.0034, 270.)
        assert absliq == 0.09822164244021624

        LiqAbsModel.model = 'R98'
        absliq = LiqAbsModel.liquid_water_absorption(0.05, 183.0034, 270.)
        assert absliq != 0.09822164244021624
        
    def test_h20continum(self):
        H2OAbsModel.model = 'R23SD'
        H2OAbsModel.set_ll()
        assert H2OAbsModel.h2oll.ctr[0] == 296.0
        
        cs = H2OAbsModel().h2o_continuum(183, 1.09, 1)
        assert_allclose(cs, np.array([3.791294e-08]))


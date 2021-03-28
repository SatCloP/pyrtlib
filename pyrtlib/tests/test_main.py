from unittest import TestCase
import os
import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
import pytest
from pathlib import Path
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import tb_cloud_rte
from pyrtlib.utils import ppmv2gkg, mr2rh

TEST_DIR = Path(__file__).parent
DATA_DIR = os.path.join(TEST_DIR, 'data')

@pytest.mark.datafiles(DATA_DIR)
def test_tb_cloud_rte_with_rose19sd(datafiles):
    z, p, d, t, md = atmp.gl_atm(atmp.TROPICAL)

    gkg = ppmv2gkg(md[:, 0], atmp.H20)
    rh = mr2rh(p, t, gkg)[0] / 100

    ang = np.array([90.])
    frq = np.arange(20, 201, 1)
    nf = len(frq)

    denliq = np.zeros(z.shape)
    denice = np.zeros(z.shape)
    cldh = np.zeros((2, 0))

    df = tb_cloud_rte(z, p, t, rh, denliq, denice, cldh, frq, ang, absmdl='rose19sd', ray_tracing=True,
                      from_sat=True)
    df = df.set_index(frq)
    df_expected = pd.read_csv(os.path.join(datafiles, "test_rose19sd_tb.csv"))
    assert_allclose(df.tbtotal, df_expected.tbtotal, atol=0)

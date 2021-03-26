# -*- coding: utf-8 -*-
"""
Main script !!!(right now only for testing)!!!.
"""

import numpy as np
import pandas as pd
from .absmodel import AbsO2Model, AbsH2OModel, AbsN2Model, AbsLiqModel
from .rte import RTEquation
from .utils import arange


def tb_cloud_rte(z=None, p=None, tk=None, rh=None, denliq=None, denice=None, cldh=None, frq=None, angles=None,
                 absmdl=None, ray_tracing=None, es=None, *args, **kwargs):

    AbsO2Model.model = absmdl
    AbsH2OModel.model = absmdl
    AbsN2Model.model = absmdl
    AbsLiqModel.model = absmdl
    RTEquation.from_sat = es

    # Settings
    nl = len(z)
    nf = len(frq)
    nang = len(angles)
    ncld = len(cldh[1, :])
    ice = 0

    # Allocation
    SPtaudry = np.zeros((nf, nang))
    SPtauwet = np.zeros((nf, nang))
    SPtauliq = np.zeros((nf, nang))
    SPtauice = np.zeros((nf, nang))
    Ptaudry = np.zeros((nf, nang, nl))
    Ptaulay = np.zeros((nf, nang, nl))
    Ptauwet = np.zeros((nf, nang, nl))
    Ptauliq = np.zeros((nf, nang, nl))
    Ptauice = np.zeros((nf, nang, nl))
    tbtotal = np.zeros((nf, nang))
    tbatm = np.zeros((nf, nang))
    tmr = np.zeros((nf, nang))
    tmrcld = np.zeros((nf, nang))
    srho = np.zeros((1, nang))
    swet = np.zeros((1, nang))
    sdry = np.zeros((1, nang))
    sliq = np.zeros((1, nang))
    sice = np.zeros((1, nang))
    # ... convert height profile to (km above antenna height) ...
    z0 = z[0]
    z = z - z0
    # ... compute vapor pressure and vapor density ...
    e, rho = RTEquation.vapor_xxx(tk, rh, ice)
    # ... convert cloud base and cloud top to (km above antenna height) ...
    # ... compute (beglev) and (endlev) ...
    cldh = cldh - z0
    # TODO: check if shape is correct
    beglev = np.zeros((1, nang))
    endlev = np.zeros((1, nang))
    for l in arange(0, ncld-1).reshape(-1):
        for i in arange(1, nl).reshape(-1):
            if z[i] == cldh[0, l]: beglev[l] = i
            if z[i] == cldh[1, l]: endlev[l] = i

    # ... compute refractivity ...
    dryn, wetn, refindx = RTEquation.refract_xxx(p, tk, e)
    for k in arange(0, nang-1).reshape(-1):
        # ... Compute distance between each level (ds) ...
        if ray_tracing:
            ds = RTEquation.ray_trac_xxx(z, refindx, angles[k], z0)
        else:
            amass = 1 / np.sin(np.dot(angles[k], np.pi) / 180)
            # TODO: check if array shape is correct
            ds = np.concatenate([[0], [np.dot(np.diff(z), amass)]])
        # ds = [0; diff(z)]; # in alternative simple diff of z
        # ... Integrate over path (ds) ...
        srho[k], _ = RTEquation.exp_int_xxx(1, rho, ds, 0, nl, 0.1)
        swet[k], _ = RTEquation.exp_int_xxx(1, wetn, ds, 0, nl, 0.1)
        sdry[k], _ = RTEquation.exp_int_xxx(1, dryn, ds, 0, nl, 0.1)
        if ncld > 0:
            sliq[k] = RTEquation.cld_int_xxx(denliq, ds, beglev, endlev)
            sice[k] = RTEquation.cld_int_xxx(denice, ds, beglev, endlev)
        # ... handle each frequency ...
        # this are based on NOAA RTE fortran routines
        for j in arange(0, nf-1).reshape(-1):
            # Rosenkranz, personal communication, 2019/02/12 (email)
            awet, adry = RTEquation.clr_abs_xxx(p, tk, e, frq[j])
            aliq, aice = RTEquation.cld_abs_xxx(tk, denliq, denice, frq[j])

            SPtauwet[j, k], Ptauwet[j, k, :] = RTEquation.exp_int_xxx(1, awet, ds, 0, nl, 1)
            SPtaudry[j, k], Ptaudry[j, k, :] = RTEquation.exp_int_xxx(1, adry, ds, 0, nl, 1)
            SPtauliq[j, k], Ptauliq[j, k, :] = RTEquation.exp_int_xxx(0, aliq, ds, 0, nl, 1)
            SPtauice[j, k], Ptauice[j, k, :] = RTEquation.exp_int_xxx(0, aice, ds, 0, nl, 1)
            Ptaulay[j, k, :] = Ptauwet[j, k, :] + Ptaudry[j, k, :] + Ptauice[j, k, :] + Ptauliq[j, k, :]
            # [boftotl,boftatm,boftmr,PSPtauprof,hvk] = Planck_xxx(frq(j),tk,Ptaulay(j,k,:));
            boftotl, boftatm, boftmr, PSPtauprof, hvk, _, _ = RTEquation.planck_xxx(frq[j], tk, Ptaulay[j, k, :])
            if ncld > 0:
                tmrcld[j, k] = RTEquation.cld_tmr_xxx(beglev[0], endlev[0], hvk, PSPtauprof, boftatm)
            # ... assign output values ...
            tbtotal[j, k] = RTEquation.bright_xxx(hvk, boftotl)
            tbatm[j, k] = RTEquation.bright_xxx(hvk, boftatm[nl-1])
            tmr[j, k] = RTEquation.bright_xxx(hvk, boftmr)

    df = pd.DataFrame({'tbtotal': tbtotal.T[0],
                       'tbatm': tbatm.T[0],
                       'tmr': tmr.T[0],
                       'tauwet': SPtauwet.T[0],
                       'taudry': SPtaudry.T[0],
                       'tauliq': SPtauliq.T[0],
                       'tauice': SPtauice.T[0]
                       })

    # 'tmrcld': tmrcld.T[0],
    # 'taulaywet': Ptauwet.T[0],
    # 'taulaydry': Ptaudry.T[0],
    # 'taulayliq': Ptauliq.T[0],
    # 'taulayice': Ptauice.T[0],
    # 'srho': srho.T[0],
    # 'swet': swet.T[0],
    # 'sdry': sdry.T[0],
    # 'sliq': sliq.T[0],
    # 'sice': sice.T[0]
    return df

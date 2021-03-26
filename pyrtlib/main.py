# -*- coding: utf-8 -*-
"""
Main script !!!(right now only for testing)!!!.
"""

import numpy as np
import pandas as pd
from .absmodel import O2AbsModel, H2OAbsModel, N2AbsModel, LiqAbsModel
from .rte import RTEquation
from .utils import arange


def tb_cloud_rte(z, p, tk, rh, denliq, denice, cldh, frq, angles, absmdl, ray_tracing=True, from_sat=True, *args, **kwargs):
    """[summary]

    Args:
        z ([type]): [description]
        p ([type]): [description]
        tk ([type]): [description]
        rh ([type]): [description]
        denliq ([type]): [description]
        denice ([type]): [description]
        cldh ([type]): [description]
        frq ([type]): [description]
        angles ([type]): [description]
        absmdl ([type]): [description]
        ray_tracing (bool, optional): [description]. Defaults to True.
        from_sat (bool, optional): [description]. Defaults to True.

    Returns:
        [type]: [description]
    """

    O2AbsModel.model = absmdl
    H2OAbsModel.model = absmdl
    N2AbsModel.model = absmdl
    LiqAbsModel.model = absmdl
    RTEquation.from_sat = from_sat

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
    e, rho = RTEquation.vapor(tk, rh, ice)
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
    dryn, wetn, refindx = RTEquation.refractivity(p, tk, e)
    for k in arange(0, nang-1).reshape(-1):
        # ... Compute distance between each level (ds) ...
        if ray_tracing:
            ds = RTEquation.ray_tracing(z, refindx, angles[k], z0)
        else:
            amass = 1 / np.sin(np.dot(angles[k], np.pi) / 180)
            # TODO: check if array shape is correct
            ds = np.concatenate([[0], [np.dot(np.diff(z), amass)]])
        # ds = [0; diff(z)]; # in alternative simple diff of z
        # ... Integrate over path (ds) ...
        srho[k], _ = RTEquation.exponential_integration(1, rho, ds, 0, nl, 0.1)
        swet[k], _ = RTEquation.exponential_integration(1, wetn, ds, 0, nl, 0.1)
        sdry[k], _ = RTEquation.exponential_integration(1, dryn, ds, 0, nl, 0.1)
        if ncld > 0:
            sliq[k] = RTEquation.cloud_integrated_density(denliq, ds, beglev, endlev)
            sice[k] = RTEquation.cloud_integrated_density(denice, ds, beglev, endlev)
        # ... handle each frequency ...
        # this are based on NOAA RTE fortran routines
        for j in arange(0, nf-1).reshape(-1):
            # Rosenkranz, personal communication, 2019/02/12 (email)
            awet, adry = RTEquation.clearsky_absorption(p, tk, e, frq[j])
            aliq, aice = RTEquation.cloudy_absorption(tk, denliq, denice, frq[j])

            SPtauwet[j, k], Ptauwet[j, k, :] = RTEquation.exponential_integration(1, awet, ds, 0, nl, 1)
            SPtaudry[j, k], Ptaudry[j, k, :] = RTEquation.exponential_integration(1, adry, ds, 0, nl, 1)
            SPtauliq[j, k], Ptauliq[j, k, :] = RTEquation.exponential_integration(0, aliq, ds, 0, nl, 1)
            SPtauice[j, k], Ptauice[j, k, :] = RTEquation.exponential_integration(0, aice, ds, 0, nl, 1)
            Ptaulay[j, k, :] = Ptauwet[j, k, :] + Ptaudry[j, k, :] + Ptauice[j, k, :] + Ptauliq[j, k, :]
            # [boftotl,boftatm,boftmr,PSPtauprof,hvk] = Planck_xxx(frq(j),tk,Ptaulay(j,k,:));
            boftotl, boftatm, boftmr, PSPtauprof, hvk, _, _ = RTEquation.planck(frq[j], tk, Ptaulay[j, k, :])
            if ncld > 0:
                tmrcld[j, k] = RTEquation.cloud_radiating_temperature(beglev[0], endlev[0], hvk, PSPtauprof, boftatm)
            # ... assign output values ...
            tbtotal[j, k] = RTEquation.bright(hvk, boftotl)
            tbatm[j, k] = RTEquation.bright(hvk, boftatm[nl - 1])
            tmr[j, k] = RTEquation.bright(hvk, boftmr)

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

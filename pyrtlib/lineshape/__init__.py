"""This version correspond to 2018/07/23 package, where h2o_list.asc is dated July 16 2018

    REFERENCES FOR MEASUREMENTS (freq, W0,D in GHZ/bar, XW, XD, A, W2) updated Feb. 20, 2019
    REFERENCES FOR MEASUREMENTS (freq, W0air, W0self, Dair, Dself in GHZ/bar, Xair, Xself, XDair, XDself, W2)

    References
    ----------
    .. [1] M. Koshelev et al., JQSRT v.205, pp. 51-58 (2018)
    .. [2] M. Koshelev (private comm., 2019)
    .. [3] G. Golubiatnikov, J. MOLEC. SPEC. vol. 230, pp.196-198 (2005)
    .. [4] M. Koshelev et al., J. Molec. Spec. v.241, pp.101-108 (2007)
    .. [5] J.-M. Colmont et al.,J. Molec. Spec. v.193, pp.233-243 (1999)
    .. [6] M. Tretyakov et al, JQSRT v.114 pp.109-121 (2013)
    .. [7] G. Golubiatnikov et al., JQSRT v.109, pp.1828-1833 (2008)
    .. [8] V. Podobedov et al., JQSRT v.87, pp. 377-385 (2004)
    .. [9] M. Koshelev, JQSRT v.112, pp.550-552 (2011)
    .. [10] M. Tretyakov, JQSRT v.328, pp.7-26 (2016)
    .. [11] V. Payne et al.,IEEE Trans. Geosci. Rem. Sens. v.46, pp.3601-3617 (2008)
    .. [12] D. Turner et al., IEEE Trans. Geosci. Rem. Sens. v.47 pp.3326-37 (2009)

    Other parameters from HITRAN2016.
    Continuum re-adjusted for new line par. Mar. 20, 2019.
"""

from ..absmod_uncertainty import uncertainty_propagation

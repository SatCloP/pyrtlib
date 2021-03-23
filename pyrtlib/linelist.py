
# This version correspond to 2018/07/23 package, where h2o_list.asc is dated July 16 2018
    
# REFERENCES FOR MEASUREMENTS (freq, W0,D in GHZ/bar, XW, XD, A, W2) updated Feb. 20, 2019
# REFERENCES FOR MEASUREMENTS (freq, W0air, W0self, Dair, Dself in GHZ/bar, Xair, Xself, XDair, XDself, W2)
# (1) M. Koshelev et al., JQSRT v.205, pp. 51-58 (2018)
# (2) M. Koshelev (private comm., 2019)
# (3) G. Golubiatnikov, J. MOLEC. SPEC. vol. 230, pp.196-198 (2005)
# (4) M. Koshelev et al., J. Molec. Spec. v.241, pp.101-108 (2007)
# (5) J.-M. Colmont et al.,J. Molec. Spec. v.193, pp.233-243 (1999)
# (6) M. Tretyakov et al, JQSRT v.114 pp.109-121 (2013)
# (7) G. Golubiatnikov et al., JQSRT v.109, pp.1828-1833 (2008)
# (8) V. Podobedov et al., JQSRT v.87, pp. 377-385 (2004)
# (9) M. Koshelev, JQSRT v.112, pp.550-552 (2011)
# (10) M. Tretyakov, JQSRT v.328, pp.7-26 (2016)
# (11) V. Payne et al.,IEEE Trans. Geosci. Rem. Sens. v.46, pp.3601-3617 (2008)
# (12) D. Turner et al., IEEE Trans. Geosci. Rem. Sens. v.47 pp.3326-37 (2009),
# Other parameters from HITRAN2016.
# Continuum re-adjusted for new line par. Mar. 20, 2019.

# 2019/03/18 - Nico: first created

#blk = -9999; # blank
blk=copy(nan)

# molecule freq,GHz  S(296K)    B     W0air XWair W0self XWself  Dair  XDair  Dself XDself  Aair Aself W2air W2self Refs.
# FL(i)     S1(i)    B2(i)   Wair  X(i)  Wself  Xs(i)   Sair  Xh(I)  Sself  Xhs(I) Aair Aself   W2   W2S (variable names in the code) Refs.
MTX=concat([[11,22.23508,1.335e-14,2.172,2.74,0.76,13.63,1.2,- 0.033,2.6,0.814,blk,blk,blk,blk,blk],[11,183.310087,2.319e-12,0.677,3.028,0.55,15.01,0.79,- 0.073,2.0,0.112,1.43,0.0,18.3,0.406,1.499],[11,321.22563,7.657e-14,6.262,2.426,0.73,10.65,0.54,- 0.143,blk,0.278,blk,blk,blk,blk,blk],[11,325.152888,2.721e-12,1.561,2.847,0.64,13.95,0.74,- 0.013,blk,1.325,blk,blk,blk,blk,blk],[11,380.197353,2.477e-11,1.062,2.868,0.54,14.4,0.89,- 0.074,blk,0.24,blk,blk,blk,blk,blk],[11,439.150807,2.137e-12,3.643,2.055,0.69,9.06,0.52,0.051,blk,0.165,blk,blk,blk,blk,blk],[11,443.018343,4.44e-13,5.116,1.819,0.7,7.96,0.5,0.14,blk,- 0.229,blk,blk,blk,blk,blk],[11,448.001085,2.588e-11,1.424,2.612,0.7,13.01,0.67,- 0.116,blk,- 0.615,blk,blk,blk,blk,blk],[11,470.888999,8.196e-13,3.645,2.169,0.73,9.7,0.65,0.061,blk,- 0.465,blk,blk,blk,blk,blk],[11,474.689092,3.268e-12,2.411,2.366,0.71,11.24,0.64,- 0.027,blk,- 0.72,blk,blk,blk,blk,blk],[11,488.490108,6.628e-13,2.89,2.616,0.75,13.58,0.72,- 0.065,blk,- 0.36,blk,blk,blk,blk,blk],[11,556.935985,1.57e-09,0.161,3.115,0.75,14.24,1.0,0.187,blk,- 1.693,blk,blk,blk,blk,blk],[11,620.700807,1.7e-11,2.423,2.468,0.79,11.94,0.75,0.0,blk,0.687,0.92,blk,blk,blk,blk],[11,658.006072,9.033e-13,7.921,3.154,0.73,13.84,1.0,0.176,blk,- 1.496,blk,blk,blk,blk,blk],[11,752.033113,1.035e-09,0.402,3.114,0.77,13.58,0.84,0.162,blk,- 0.878,blk,blk,blk,blk,blk],[11,916.171582,4.275e-11,1.461,2.695,0.79,13.55,0.48,0.0,blk,0.521,0.47,blk,blk,blk,blk]])

# continuum terms
CTR=concat([300.0,5.929e-10,3.0,1.42e-08,7.5])

# Below is from abh2o_sd.f
# Read line parameters; units: GHz, Hz*cm^2, MHz/mb
REFTLINE=296.0

FL=MTX(arange(),2)

S1=MTX(arange(),3)

B2=MTX(arange(),4)

W0=MTX(arange(),5) / 1000.0

X=MTX(arange(),6)

W0S=MTX(arange(),7) / 1000.0

XS=MTX(arange(),8)

SH=MTX(arange(),9) / 1000.0

XH=MTX(arange(),10)

SHS=MTX(arange(),11) / 1000.0

XHS=MTX(arange(),12)

Aair=MTX(arange(),13)

Aself=MTX(arange(),14)

W2=MTX(arange(),15) / 1000.0

W2S=MTX(arange(),16) / 1000.0

# Replace non-existing shifting parameters with broadening parameters
#indx = find(XH==blk); XH(indx) = X(indx);
#indx = find(XHS==blk); XHS(indx) = XS(indx);
indx=find(isnan(XH))
XH[indx]=X(indx)
indx=find(isnan(XHS))
XHS[indx]=XS(indx)
# Replace non-existing Aair Aself parameters with zero (to be agreed with Phil)
#indx = find(Aair==blk); Aair(indx) = 0;
#indx = find(Aself==blk); Aself(indx) = 0;
indx=find(isnan(Aair))
Aair[indx]=0
indx=find(isnan(Aself))
Aself[indx]=0
# Also for W2 and W2S (to be agreed with Phil)
#indx = find(W2==blk); W2(indx) = 0;
#indx = find(W2S==blk); W2S(indx) = 0;
indx=find(isnan(W2))
W2[indx]=0
indx=find(isnan(W2S))
W2S[indx]=0
# Read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON=CTR(1)
CF=CTR(2)
XCF=CTR(3)
CS=CTR(4)
XCS=CTR(5)

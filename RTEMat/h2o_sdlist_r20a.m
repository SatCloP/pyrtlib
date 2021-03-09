% From /Users/Nico/SFTWR/MPM/PWR2020/20200625_revised183/h2o_sdlist.asc
%
% Two modifications wrt h2o_sdlist_r20.m
% 1) values consistent with manuscript by Koshelev et al. 2020
% 2) added delta2 values (hard-coded in Phil's routines)
%
% REFERENCES FOR MEASUREMENTS (freq, W0,D in GHZ/bar, XW, XD, A, W2) updated June. 25, 2020.
% (1) M. Koshelev et al., JQSRT v.205, pp. 51-58 (2018)
% (2) V. Payne et al.,IEEE Trans. Geosci. Rem. Sens. v.46, pp.3601-3617 (2008)
% (3) G. Golubiatnikov, J. MOLEC. SPEC. vol. 230, pp.196-198 (2005)
% (4) M. Koshelev et al., J. Molec. Spec. v.241, pp.101-108 (2007)
% (5) J.-M. Colmont et al.,J. Molec. Spec. v.193, pp.233-243 (1999)
% (6) M. Tretyakov et al, JQSRT v.114 pp.109-121 (2013)
% (7) G. Golubiatnikov et al., JQSRT v.109, pp.1828-1833 (2008)
% (8) V. Podobedov et al., JQSRT v.87, pp. 377-385 (2004)
% (9) M. Koshelev, JQSRT v.112, pp.550-552 (2011)
% (10) M. Tretyakov, JQSRT v.328, pp.7-26 (2016)
% (11) D. Turner et al., IEEE Trans. Geosci. Rem. Sens. v.47 pp.3326-37 (2009), re-adjusted for new line par. June 25, 2020.
% (12) M. Koshelev et al. (manuscript, 2020)
% Other parameters from HITRAN2016.
%
%
% 2020/06/29 - Nico: first created 

%blk = -9999; % blank
blk = nan; % blank

% molecule freq,GHz  S(296K)    B     W0air XWair W0self XWself  Dair  XDair  Dself XDself  Aair Aself W2air W2self Refs.
%          FL(i)     S1(i)    B2(i)   Wair  X(i)  Wself  Xs(i)   Sair  Xh(I)  Sself  Xhs(I) Aair Aself   W2   W2S (variable names in the code) Refs.
MTX = [...
  11   22.235080  0.1335E-13  2.172  2.740  0.76  13.63  1.20  -0.033   2.6   0.814   blk   blk   blk 0.435 1.91   ; ...% 1,2,10
  11  183.310087  0.2319E-11  0.677  3.033  0.62  15.01  0.82  -.074    1.8   0.136  0.98   0.   12.6 0.413 1.48  ; ...% 12
  11  321.225630  0.7657E-13  6.262  2.426  0.73  10.65  0.54  -0.143   blk   0.278   blk   blk   blk   blk   blk  ; ...% 4
  11  325.152888  0.2721E-11  1.561  2.847  0.64  13.95  0.74  -0.013   blk   1.325   blk   blk   blk   blk   blk  ; ...% 4,5
  11  380.197353  0.2477E-10  1.062  2.868  0.54  14.40  0.89  -0.074   blk   0.240   blk   blk   blk   blk   blk  ; ...% 6
  11  439.150807  0.2137E-11  3.643  2.055  0.69   9.06  0.52   0.051   blk   0.165   blk   blk   blk   blk   blk  ; ...% 6
  11  443.018343  0.4440E-12  5.116  1.819  0.70   7.96  0.50   0.140   blk  -0.229   blk   blk   blk   blk   blk  ; ...% 6
  11  448.001085  0.2588E-10  1.424  2.612  0.70  13.01  0.67  -0.116   blk  -0.615   blk   blk   blk   blk   blk  ; ...% 6
  11  470.888999  0.8196E-12  3.645  2.169  0.73   9.70  0.65   0.061   blk  -0.465   blk   blk   blk   blk   blk  ; ...% 6
  11  474.689092  0.3268E-11  2.411  2.366  0.71  11.24  0.64  -0.027   blk  -0.720   blk   blk   blk   blk   blk  ; ...% 6
  11  488.490108  0.6628E-12  2.890  2.616  0.75  13.58  0.72  -0.065   blk  -0.360   blk   blk   blk   blk   blk  ; ...% 6
  11  556.935985  0.1570E-08  0.161  3.115  0.75  14.24  1.00   0.187   blk  -1.693   blk   blk   blk   blk   blk  ; ...% 7
  11  620.700807  0.1700E-10  2.423  2.468  0.79  11.94  0.75     0.0   blk   0.687  0.92   blk   blk   blk   blk  ; ...% 8 blanks in Dair are set to 0 (confirmed by Phil)
  11  658.006072  0.9033E-12  7.921  3.154  0.73  13.84  1.00   0.176   blk  -1.496   blk   blk   blk   blk   blk  ; ...% 7 
  11  752.033113  0.1035E-08  0.402  3.114  0.77  13.58  0.84   0.162   blk  -0.878   blk   blk   blk   blk   blk  ; ...% 9
  11  916.171582  0.4275E-10  1.461  2.695  0.79  13.55  0.48     0.0   blk   0.521  0.47   blk   blk   blk   blk  ];   % 8 blanks in Dair are set to 0 (confirmed by Phil)

% continuum terms
CTR = [300.          5.919E-10  3.        1.42E-8    7.5];                                                              % 11

% Below is from abh2o_sd.f
% Read line parameters; units: GHz, Hz*cm^2, MHz/mb
REFTLINE = 296.0; % reference T for lines
FL     = MTX(:,2);     % LINE FREQUENCIES [GHz]
S1     = MTX(:,3);     % LINE INTENSITIES AT 296K [Hz*cm^2]
B2     = MTX(:,4);     % T COEFF. OF INTENSITIES [unitless]
W0     = MTX(:,5)/1e3; % AIR-BROADENED WIDTH PARAMETERS AT REFTLINE [MHz/mb -> GHz/mb]
X      = MTX(:,6);     % T-EXPONENT OF AIR-BROADENING [unitless]
W0S    = MTX(:,7)/1e3; % SELF-BROADENED WIDTH PARAMETERS AT REFTLINE [MHz/mb -> GHz/mb]
XS     = MTX(:,8);     % T-EXPONENT OF SELF-BROADENING [unitless]
SH     = MTX(:,9)/1e3; % AIR-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
XH     = MTX(:,10);    % T-EXPONENT OF AIR-SHIFTING [unitless]
SHS    = MTX(:,11)/1e3;% SELF-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
XHS    = MTX(:,12);    % T-EXPONENT OF SELF-SHIFTING [unitless]
Aair   = MTX(:,13);    % 
Aself  = MTX(:,14);    % 
W2     = MTX(:,15)/1e3;% 
W2S    = MTX(:,16)/1e3;% 

% Replace non-existing shifting parameters with broadening parameters
%indx = find(XH==blk); XH(indx) = X(indx);
%indx = find(XHS==blk); XHS(indx) = XS(indx);
indx = find(isnan(XH)); XH(indx) = X(indx);
indx = find(isnan(XHS)); XHS(indx) = XS(indx);

% Replace non-existing Aair Aself parameters with zero (to be agreed with Phil)
%indx = find(Aair==blk); Aair(indx) = 0;
%indx = find(Aself==blk); Aself(indx) = 0;
indx = find(isnan(Aair)); Aair(indx) = 0;
indx = find(isnan(Aself)); Aself(indx) = 0;
% Also for W2 and W2S (to be agreed with Phil)
%indx = find(W2==blk); W2(indx) = 0;
%indx = find(W2S==blk); W2S(indx) = 0;
indx = find(isnan(W2)); W2(indx) = 0;
indx = find(isnan(W2S)); W2S(indx) = 0;

% Read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON = CTR(1);
CF      = CTR(2);
XCF     = CTR(3);
CS      = CTR(4);
XCS     = CTR(5);

% Values for SD delta2 (from Koshelev et al., 2020) - hard-coded in Phil's routines, e.g.:
% /Users/Nico/SFTWR/MPM/PWR2020/20200625_revised183/abh2o_sd+d2.f
% GHz/mb, for line 2 (183 GHz)
D2AIR = -.16E-4; % -0.021 MHz/torr (Koshelev et al., 2020) * mb2torr (0.750062) * MHz2GHz (1e-3) = -1.575e-05 GHz/mb
D2SELF = .16E-3; % 0.21 MHz/torr (Koshelev et al., 2020) * mb2torr (0.750062) * MHz2GHz (1e-3) = 1.575e-04 GHz/mb

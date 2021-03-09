%       REAL FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
% c     name- o2abs_19.f  language- Fortran 77
% C
% C     RETURNS POWER ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
% C     IN NEPERS/KM.  MULTIPLY O2ABS2 BY 4.343 TO CONVERT TO DB/KM.
% C
% C      5/1/95  P. Rosenkranz 
% C      11/5/97  P. Rosenkranz - 1- line modification.
% c      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
% c      8/21/02  pwr - revised width at 425
% c      3/20/03  pwr - 1- line mixing and width revised
% c      9/29/04  pwr - new widths and mixing, using HITRAN intensities
% c                     for all lines
% c      6/12/06  pwr - chg. T dependence of 1- line to 0.8
% C      10/14/08 pwr - moved isotope abundance back into intensities, 
% c                     added selected O16O18 lines.
% c      5/30/09  pwr - remove common block, add weak lines.
% c      12/18/14 pwr - adjust line broadening due to water vapor.
% C      9/29/18  pwr - 2nd-order line mixing
% c      8/20/19  pwr - adjust intensities according to Koshelev meas.
% C
%       IMPLICIT NONE
% C
% C     ARGUMENTS:
%       REAL TEMP,PRES,VAPDEN,FREQ
%       INTEGER MIXING
% C
% C     NAME    UNITS    DESCRIPTION        VALID RANGE
% C
% C     TEMP    KELVIN   TEMPERATURE        245 to 334 (lab range)
% C     PRES   MILLIBARS PRESSURE           3 TO 1000
% C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
% C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
% C     FREQ    GHZ      FREQUENCY          0 TO 900
% C
% C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
% C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
% C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993)
% C       (http://hdl.handle.net/1721.1/68611).
% c     G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, 
% c      pp.282-287 (2003).
% C     M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
% C     M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
% c     B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
% c     D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
% c     D.S. Makarov et al, JQSRT v.112, pp.1420-1428 (2011). 
% c     M.A. Koshelev et al, JQSRT, v.154, pp.24-27 (2015).
% c     M.A. Koshelev et al, JQSRT, v.169, pp.91-95 (2016).
% c     M.A. Koshelev et al, JQSRT, v.196, pp.78?86 (2017).
% c     D.S. Makarov et al, JQSRT v.243 (March 2020) doi:10.1016/j.jqsrt.2019.106798
% C     Line intensities from HITRAN2004.
% c     Non-resonant intensity from JPL catalog.
% 
% c     notes:
% c     1. The mm line-width coefficients are from Tretyakov et al 2005,
% c        Makarov et al 2008, and Koshelev et al 2016;
% c        submm line-widths from Golubiatnikov & Krupnov, except 
% c        234-GHz line width from Drouin.
% c        Mixing coeff. from Makarov's 2018 revision.
% c     2. The same temperature dependence (X) is used for submillimeter 
% c        line widths as in the 60 GHz band: (1/T)**X (Koshelev et al 2016)
% c     3. The sign of DNU in the shape factor is corrected.
% C
%
% Nico 2020/03/06 *********************************************************
% 
% From: /Users/Nico/SFTWR/MPM/PWR2020/absorption/o2abs_19.f

function [npp] = o2abs_rosen19_xxx(pdrykpa,vx,ekpa,frq);

% NB: NL=49
%C      LINES ARE ARRANGED 1-,1+,...33-,33+ IN SPIN-ROTATION SPECTRUM;
%c      BY FREQUENCY IN SUBMM SPECTRUM.
%Nico F[GHz]
      F = [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,...
            59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,...
            56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,...
            55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,...
            53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,...
            52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310,...
            50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630,...
            487.2493, 566.8956, 715.3929, 731.1866,...
            773.8395, 834.1455, 895.0710];

%Nico S(T_0)[Hz*cm2]
      S300 = [0.2906E-14,0.7957E-15,0.2444E-14,0.2194E-14,...
              0.3301E-14,0.3243E-14,0.3664E-14,0.3834E-14,...
              0.3588E-14,0.3947E-14,0.3179E-14,0.3661E-14,...
              0.2590E-14,0.3111E-14,0.1954E-14,0.2443E-14,...
              0.1373E-14,0.1784E-14,0.9013E-15,0.1217E-14,...
              0.5545E-15,0.7766E-15,0.3201E-15,0.4651E-15,...
              0.1738E-15,0.2619E-15,0.8880E-16,0.1387E-15,...
              0.4272E-16,0.6923E-16,0.1939E-16,0.3255E-16,...
              0.8301E-17,0.1445E-16,0.3356E-17,0.6049E-17,...
              0.1280E-17,0.2394E-17,...
              0.3287E-16,0.6463E-15,0.1334E-16,0.7049E-14,...
              0.3011E-14,0.1797E-16,0.1826E-14,0.2193E-16,...
              0.1153E-13,0.3974E-14,0.2512E-16];
 
%Nico (Elow+hf)/kT_0 [unitless]
      BE = [.010, .014, .083, .083, .207, .207, .387, .387, .621, .621,...
            .910, .910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.109, 2.618, 2.618,...
           3.182, 3.182, 3.800, 3.800, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819, 6.819,... 
           7.709, 7.709, 8.653, 8.653, 9.651, 9.651,...
            .019, .048, .045, .044, .049, .084, .145, .136, .141, .145, .201];

      % C      WIDTHS IN GHz/bar
%Nico gamma(T_0) [MHZ/mb == GHz/bar]
      WB300 = .56; %Nico: non-resonannt abs width
      X = .754;      
      W300 = [1.685, 1.703, 1.513, 1.495, 1.433, 1.408,...
              1.353, 1.353, 1.303, 1.319, 1.262, 1.265, 1.238, 1.217,...
              1.207, 1.207, 1.137, 1.137, 1.101, 1.101, 1.037, 1.038,...
              0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811, 0.764, 0.764,...
              0.717, 0.717, 0.669, 0.669,...
              1.65, 1.64, 1.64, 1.64, 1.60, 1.60, 1.60, 1.60, 1.62, 1.47, 1.47];
        

      %c 1st-order mixing coeff. in 1/bar
      Y0 = [-0.041, 0.277, -0.373, 0.560, -0.573, 0.618,...
            -0.366, 0.278, -0.089, -0.021, 0.0599, -0.152,...
             0.216, -0.293, 0.374, -0.436, 0.491, -0.542,...
             0.571, -0.613, 0.636, -0.670, 0.690, -0.718,...
             0.740, -0.763, 0.788, -0.807, 0.834, -0.849,...
             0.876, -0.887, 0.915, -0.922, 0.950, -0.955,...
             0.987, -0.988,  repmat(0.0,1,11)];
      Y1 = [0.,     0.11, -0.009,  0.007, 0.049, -0.1,...
            0.260, -0.346, 0.364, -0.422, 0.315, -0.341,...
            0.483, -0.503, 0.598, -0.610, 0.630, -0.633,...
            0.613, -0.611, 0.570, -0.564, 0.58,  -0.57,...
            0.61,  -0.60,  0.64,  -0.62,  0.65,  -0.64,...
            0.66,  -0.64,  0.66,  -0.64,  0.66,  -0.64,...
            0.65,  -0.63,  repmat(0.0,1,11)];
          
      %c 2nd-order mixing coeff. in 1/bar^2
      G0 = [ -0.000695, -0.090, -0.103, -0.239, -0.172, -0.171,...
            0.028, 0.150, 0.132, 0.170, 0.087, 0.069,...
             0.083, 0.068, 0.007, 0.016, -0.021, -0.066,...
            -0.095, -0.116, -0.118, -0.140, -0.173, -0.186,...
            -0.217, -0.227, -0.234, -0.242, -0.266, -0.272,...
            -0.301, -0.304, -0.334, -0.333, -0.362, -0.358,...
            -0.348, -0.344,  repmat(0.0,1,11)];
      G1 = [ 0., -0.042, 0.004, 0.025, 0.083, 0.167, ...
             0.178, 0.223, 0.054, 0.003, 0.002, -0.044,...
             -0.019, -0.054, -0.177, -0.208, -0.294, -0.334,...
            -0.368, -0.386, -0.374, -0.384, -0.387, -0.389,...
            -0.423, -0.422, -0.46, -0.46, -0.51, -0.50,...
            -0.55, -0.53, -0.58, -0.56, -0.62, -0.59,...
            -0.68, -0.65, repmat(0.0,1,11)];

      % c dnu in GHz/bar^2  
      DNU0 = [ -0.00028, 0.00596, -0.01950, 0.032, -0.0475, 0.0541,...
              -0.0232, 0.0155, 0.0007, -0.0086, -0.0026, -0.0013,...
              -0.0004, -0.002, 0.005, -0.007, 0.007, -0.008,...
               0.006, -0.007, 0.006, -0.006, 0.005, -0.0049,...
              0.0040, -0.0041, 0.0036, -0.0037, 0.0033, -0.0034,...
               0.0032, -0.0032, 0.0030, -0.0030, 0.0028, -0.0029,...
               0.0029, -0.0029, repmat(0.0,1,11)];
       DNU1 = [ -0.00037, 0.0086, -0.013, 0.019, -0.026, 0.027,...
              0.005, -0.014, 0.012, -0.018, -0.015, 0.015,...
               0.003, -0.004, 0.012, -0.013, 0.012, -0.012,...
               0.009, -0.009, 0.002, -0.002, 0.0005, -0.0005,...
              0.002, -0.002, 0.002, -0.002, 0.002, -0.002,...
              0.002, -0.002, 0.002, -0.002, 0.001, -0.001,...
              0.0004, -0.0004,  repmat(0.0,1,11)];

% CYH*** add the following lines *************************
      db2np = log(10.) * 0.1;
      rvap = 0.01 * 8.314510 / 18.01528;
      factor = .182 * frq;
      TEMP = 300./vx;
      PRES = (pdrykpa+ekpa)*10.;
      VAPDEN = ekpa*10./(rvap*TEMP);
      FREQ = frq;
% CYH*****************************************************

      TH = 300./TEMP;
      TH1 = TH-1.;
      B = TH^X;
      PRESWV = VAPDEN*TEMP/216.68;            %Nico hPa
      PRESDA = PRES -PRESWV;                  %Nico hPa
      DEN = .001*(PRESDA*B + 1.2*PRESWV*TH);  %Nico bar
      DFNR = WB300*DEN;                       %Nico GHz/bar * bar = GHz
      PE2 = DEN*DEN;

%Nico intensities of the non-resonant transitions for O16-O16 and O16-O18, from JPL's line compilation
%c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
      SUM = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));

      NLINES = length(F);
      for K = 1:NLINES
         Y = DEN * ( Y0(K) + Y1(K)*TH1 );
         DNU = PE2 * (DNU0(K) + DNU1(K)*TH1 );
         GFAC = 1. + PE2 * ( G0(K) + G1(K)*TH1 );
         DF = W300(K)*DEN;
         STR = S300(K)*exp(-BE(K)*TH1);
         DEL1 = FREQ - F(K) - DNU;
         DEL2 = FREQ + F(K) + DNU;
         D1 = DEL1*DEL1 + DF*DF;
         D2 = DEL2*DEL2 + DF*DF;
         SF1 = (DF*GFAC + DEL1*Y)/D1;
         SF2 = (DF*GFAC - DEL2*Y)/D2;
         SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))^2;
      end

%c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11
      O2ABS = 1.6097E11*SUM*PRESDA*TH^3;
      O2ABS = 1.004 * max(O2ABS,0.); % increase absorption to match Koshelev2017
      
% CYH *** ********************************************************
% C   separate the equ. into line and continuum
% C   terms, and change the units from np/km to ppm

      npp = O2ABS; % this includes the non-resonant term below (ncpp is no longer added, but it's computed just for reference)

%Nico intensities of the non-resonant transitions for O16-O16 and O16-O18, from JPL's line compilation
%c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
%      ncpp = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));
%c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11      
%Nico a/(pi*k*T_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont for units conversion (Pa->hPa, Hz->GHz)
%Nico Pa2hPa=1e-2; Hz2GHz=1e-9; m2cm=1e2; m2km=1e-3; Pa2hPa^-1 * Hz2GHz * m2cm^-2 * m2km^-1 = 1e-8
%Nico TH^3 = TH(from ideal gas law 2.13) * TH(from the MW approx of stimulated emission 2.16 vs. 2.14) * TH(from the partition sum 2.20)
%      ncpp = 1.6097E11*ncpp*PRESDA*TH^3; %Nico: n/pi*SUM0
      
% C     change the units from np/km to ppm
      npp = (npp /db2np)/factor; % this includes ncpp; do NOT sum ncpp
%      ncpp = (ncpp / db2np)/factor;

% CYH ************************************************************
         
return

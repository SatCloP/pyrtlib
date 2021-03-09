% C************************************
% C The interface of the following original fuction is changed to match
% C the interface of the ETL routines.    Yong Han, 2000.
% CYH*****************************************************************
%       REAL FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
% C  Copyright (c) 2010 Massachusetts Institute of Technology
% C
% C     RETURNS POWER ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
% C     IN NEPERS/KM.  MULTIPLY O2ABS BY 4.343 TO CONVERT TO DB/KM.
% C
% C      12/14/10  P. Rosenkranz 
% C
%       IMPLICIT NONE
% C
% C     ARGUMENTS:
%       REAL TEMP,PRES,VAPDEN,FREQ
% C
% C     NAME    UNITS    DESCRIPTION        VALID RANGE
% C
% C     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
% c                                          valid for atmosphere
% C     PRES   MILLIBARS PRESSURE           3 TO 1000
% C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
% C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
% C     FREQ    GHZ      FREQUENCY          0 TO 900
% C
% C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
% C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
% C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
% c     G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, 
% c      pp.282-287 (2003).
% C     M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
% C     M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
% c     A. I. Meshkov & F.C. De Lucia, JQSRT v.108, pp. 256-276 (2007).
% c     D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
% c     D.S. Makarov et al, JQSRT v.112, pp.1420-1428 (2011). 
% C     Line intensities from HITRAN2004.
% c     Non-resonant intensity from JPL catalog.
% c     mixing coeff. from Makarov et al 2011.
% c     note:
% c     1. The mm line-width coefficients are from Tretyakov et al 2005
% c        and Makarov et al 2008;
% c        submm line-widths from Golubiatnikov & Krupnov, except 
% c        234-GHz line width from Meshkov & De Lucia.
% c     2. The same temperature dependence (X) is used for submillimeter 
% c        line widths as in the 60 GHz band: (1/T)**X
% c     3. The sign of DNU in the shape factor is corrected.

%      subroutine O2N2_xxx (pdrykpa,vx,ekpa,frq,npp,ncpp)
function [npp,ncpp] = o2n2_rosen16_xxx(pdrykpa,vx,ekpa,frq);

% Nico 2016/11/30 *********************************************************
% Here I imported the code I got from P. Rosenkranz on 2016/08/10, adapting
% to RTE. 

% NB: NL=49
%C      LINES ARE ARRANGED 1-,1+,...33-,33+ IN SPIN-ROTATION SPECTRUM;
%c      BY FREQUENCY IN SUBMM SPECTRUM.
      F = [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,...
            59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,...
            56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,...
            55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,...
            53.5958, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,...
            52.0214, 67.3696, 51.5034, 67.9009, 50.9877, 68.4310,...
            50.4742, 68.9603, 233.9461, 368.4982, 401.7398, 424.7630,...
            487.2493, 566.8956, 715.3929, 731.1866,...
            773.8395, 834.1455, 895.0710];

     
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

          
      BE = [.010, .014, .083, .083, .207, .207, .387, .387, .621, .621,...
            .910, .910, 1.255, 1.255, 1.654, 1.654, 2.109, 2.109, 2.618, 2.618,...
           3.182, 3.182, 3.800, 3.800, 4.474, 4.474, 5.201, 5.201, 5.983, 5.983, 6.819, 6.819,... 
           7.709, 7.709, 8.653, 8.653, 9.651, 9.651,...
            .019, .048, .045, .044, .049, .084, .145, .136, .141, .145, .201];

      % C      WIDTHS IN MHZ/MB
      WB300 = .56; 
      X = .8;      
      W300 = [1.688, 1.703, 1.513, 1.491, 1.415, 1.408,...
              1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217,...
              1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038,...
              0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811, 0.764, 0.764,...
              0.717, 0.717, 0.669, 0.669,...
              2.78, 1.64, 1.64, 1.64, 1.60, 1.60, 1.60, 1.60, 1.62, 1.47, 1.47];
           

      %c 1st-order mixing coeff. 
      Y0 = [-0.041,  0.265, -0.354,   0.53,  -0.547,   0.589,...
            -0.359,  0.276, -0.0942, -0.0124, 0.0647, -0.155,...
             0.224, -0.299,  0.365,  -0.426,  0.466,  -0.515,...
             0.546, -0.587,  0.608,  -0.642,  0.64,   -0.668,...
             0.655, -0.679,  0.669,  -0.689,  0.689,  -0.706,...
             0.715, -0.729,  0.743,  -0.754,  0.771,  -0.78,...
             0.802, -0.807,  repmat(0.0,1,11)];
      Y1 = [0.000, -0.0132, 0.0582, -0.0981,...
            0.114, -0.153, -0.0494,  0.0476, 0.186, -0.235,...
            0.611, -0.674,  0.628,  -0.662,  0.583, -0.598,...
            0.596, -0.6,    0.617,  -0.615,  0.515, -0.509,...
            0.406, -0.399,  0.348,  -0.34,   0.335, -0.327,...
            0.342, -0.332,  0.346,  -0.336,  0.344, -0.334,...
            0.338, -0.327,  0.32,   -0.309,  repmat(0.0,1,11)];
          
      %c 2nd-order mixing coeff. 
      G0 = [ 0.000,  -0.0835, -0.0947, -0.218,...
            -0.16,   -0.162,   0.0197,  0.132,   0.128,   0.166,...
             0.0883,  0.0713,  0.0735,  0.0582, -0.00369, 0.00439,...
            -0.0213, -0.0634, -0.0868, -0.105,  -0.114,  -0.134,...
            -0.163,  -0.174,  -0.186,  -0.195,  -0.179,  -0.186,...
            -0.184,  -0.189,  -0.197,  -0.2,    -0.212,  -0.214,...
            -0.227,  -0.228,  -0.215,  -0.216,  repmat(0.0,1,11)];
      G1 = [ 0.000,  0.00843, 0.0234, 0.0596, ...
             0.0326, 0.041,   0.0174, 0.039,  0.262,  0.305,...
             0.126,  0.0799, -0.156, -0.216, -0.253, -0.296,...
            -0.287, -0.324,  -0.374, -0.394, -0.393, -0.4,...
            -0.329, -0.328,  -0.261, -0.258, -0.237, -0.232,...
            -0.244, -0.239,  -0.262, -0.254, -0.278, -0.268,...
            -0.295, -0.283,  -0.354, -0.336, repmat(0.0,1,11)];
      DNU0 = [ 0.000,    0.00545, -0.0175,   0.0286,...
              -0.043,    0.0491,  -0.0229,   0.0162,  -0.000866,...
              -0.00609, -0.00115, -0.00256,  0.00221, -0.00447,...
               0.00605, -0.00746,  0.00675, -0.00762,  0.00732,...
              -0.00784,  0.00802, -0.00833,  0.00696, -0.00714,...
               0.00517, -0.00526,  0.00368, -0.00376,  0.00277,...
              -0.00283,  0.0023,  -0.00235,  0.00209, -0.00211,...
               0.00194, -0.00197,  0.00199, -0.00199, repmat(0.0,1,11)];
      DNU1 = [ 0.000,     0.0017,   -0.0013,    0.00107,...
              -0.00444,   0.00493,  -0.0529,    0.061,    -0.0193,...
               0.0168,    0.0393,   -0.0448,    0.0247,   -0.0265,...
               0.0102,   -0.0106,    0.0104,   -0.0105,    0.0163,...
              -0.0163,    0.00835,  -0.00815,   0.000922, -0.00076,...
              -0.00174,   0.00184,  -0.0015,    0.00157,  -0.000561,...
               0.000641, -0.0000939, 0.000173, -0.0000366, 0.000099,...
              -0.000156,  0.000224, -0.000708,  0.000743,  repmat(0.0,1,11)];
        
% Nico 2016/11/30 *********************************************************

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
      PRESWV = VAPDEN*TEMP/217.;
      PRESDA = PRES -PRESWV;
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH);
      DFNR = WB300*DEN;
      PE2 = DEN*DEN;

% CYH *** The following line is taken out, but **********************
% C       it is included in the ncpp term 
%c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
%      SUM = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      SUM = 0.0;

% CYH **************************************************************

      NLINES = length(F);
      for K = 1:NLINES
         DF = W300(K)*DEN;
         %Y = .001*PRES*B*(Y300(K)+V(K)*TH1);
         Y = DEN * ( Y0(K) + Y1(K)*TH1 );
         DNU = PE2 * (DNU0(K) + DNU1(K)*TH1 );
         STR = S300(K)*exp(-BE(K)*TH1);
         DEL1 = FREQ - F(K) - DNU;
         DEL2 = FREQ + F(K) + DNU;
         GFAC = 1. + PE2 * ( G0(K) + G1(K)*TH1 );
         D1 = DEL1*DEL1 + DF*DF;
         D2 = DEL2*DEL2 + DF*DF;
         SF1 = (DF*GFAC + DEL1*Y)/D1;
         SF2 = (DF*GFAC - DEL2*Y)/D2;
         SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))^2;
      end

%      O2ABS = .5034E12*SUM*PRESDA*TH^3/3.14159;
%c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11      
      O2ABS = 1.6097E11*SUM*PRESDA*TH^3;
      O2ABS = max(O2ABS,0.);

% CYH *** ********************************************************
% C   separate the equ. into line and continuum
% C   terms, and change the units from np/km to ppm

      npp = O2ABS;
%c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
      ncpp = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));
%c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11      
      ncpp = 1.6097E11*ncpp*PRESDA*TH^3;
% C    add N2 term
      ncpp = ncpp + ABSN2_ros16(TEMP,PRES,FREQ);
% C     change the units from np/km to ppm
      npp = (npp /db2np)/factor;
      ncpp = (ncpp / db2np)/factor;

% CYH ************************************************************

return
      
end
      
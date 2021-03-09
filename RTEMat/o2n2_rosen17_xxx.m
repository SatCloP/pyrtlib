% C************************************
% C The interface of the following original fuction is changed to match
% C the interface of the ETL routines.    Yong Han, 2000.
% CYH*****************************************************************
%
%       REAL FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
% C  Copyright (c) 2009 Massachusetts Institute of Technology
% C
% C     RETURNS POWER ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
% C     IN NEPERS/KM.  MULTIPLY O2ABS BY 4.343 TO CONVERT TO DB/KM.
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
% C     P.W. Rosenkranz, CHAP. 2 in ATMOSPHERIC REMOTE SENSING
% C       BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993) 
% C       (http://hdl.handle.net/1721.1/68611).
% c     G.Yu. Golubiatnikov & A.F. Krupnov, J. Mol. Spect. v.217, 
% c       pp.282-287 (2003).
% C     M.Yu. Tretyakov et al, J. Mol. Spect. v.223, pp.31-38 (2004).
% C     M.Yu. Tretyakov et al, J. Mol. Spect. v.231, pp.1-14 (2005).
% c     B.J. Drouin, JQSRT v.105, pp.450-458 (2007).
% c     D.S. Makarov et al, J. Mol. Spect. v.252, pp.242-243 (2008).
% c     M.A. Koshelev et al, JQSRT, in press (2015).
% C     line intensities from HITRAN2004.
% c     non-resonant intensity from JPL catalog.
% c
% c     note:
% c     1. The mm line-width and mixing coefficients are from Tretyakov et al;
% c        submm line-widths from Golubiatnikov & Krupnov (except 
% c        234 GHz from Drouin)
% c     2. The same temperature dependence (X) is used for submillimeter 
% c        line widths as in the 60 GHz band: (1/T)**X 
%
% Nico 2018/02/18 *********************************************************
% References:
%
% Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017

%      subroutine O2N2_xxx (pdrykpa,vx,ekpa,frq,npp,ncpp)
function [npp,ncpp] = o2n2_rosen17_xxx(pdrykpa,vx,ekpa,frq);

% Nico 2016/11/30 *********************************************************
% Here I imported the code I got from P. Rosenkranz on 2016/08/10, adapting
% to RTE. 

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

      % C      WIDTHS IN MHZ/MB
%Nico gamma(T_0) [MHZ/mb == GHz/bar]
      WB300 = .56; %Nico: non-resonannt abs width
      X = .8;      
      W300 = [1.688, 1.703, 1.513, 1.491, 1.415, 1.408,...
              1.353, 1.339, 1.295, 1.292, 1.262, 1.263, 1.223, 1.217,...
              1.189, 1.174, 1.134, 1.134, 1.089, 1.088, 1.037, 1.038,...
              0.996, 0.996, 0.955, 0.955, 0.906, 0.906, 0.858, 0.858, 0.811, 0.811, 0.764, 0.764,...
              0.717, 0.717, 0.669, 0.669,...
              1.65, 1.64, 1.64, 1.64, 1.60, 1.60, 1.60, 1.60, 1.62, 1.47, 1.47];

%Nico Y(T_0) [unitless]
      Y300 = [-0.0360, 0.2547, -0.3655,  0.5495,...
              -0.5696,  0.6181, -0.4252,  0.3517, -0.1496,  0.0430,...
               0.0640, -0.1605,  0.2906, -0.3730,  0.4169, -0.4819,...
               0.4963, -0.5481,  0.5512, -0.5931,  0.6212, -0.6558,...
               0.6920, -0.7208,  0.7312, -0.7550,  0.7555, -0.7751,...
               0.7914, -0.8073,  0.8307, -0.8431,  0.8676, -0.8761,...
               0.9046, -0.9092,  0.9416, -0.9423,...
               repmat(0.0,1,11)];

%Nico v(T_0) [unitless]
      V = [  0.0079, -0.0978,  0.0844, -0.1273,...
             0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,...
             0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,...
             0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,...
             0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,...
             0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,...
             0.6800, -0.660,   0.685,  -0.665,...
               repmat(0.0,1,11)];
             
% Nico 2016/11/30 *********************************************************

% CYH*** add the following lines *************************
      db2np = log(10.) * 0.1;
      rvap = 0.01 * 8.314510 / 18.01528;
      factor = .182 * frq;
      TEMP = 300./vx;
      PRES = (pdrykpa+ekpa)*10.; %Nico kPa to hPa
      VAPDEN = ekpa*10./(rvap*TEMP);
      FREQ = frq;
% CYH*****************************************************

      TH = 300./TEMP;
      TH1 = TH-1.;
      B = TH^X;
      PRESWV = VAPDEN*TEMP/217.;             %Nico hPa
      PRESDA = PRES -PRESWV;                 %Nico hPa
      DEN = .001*(PRESDA*B + 1.2*PRESWV*TH); %Nico bar
      DFNR = WB300*DEN;                      %Nico GHz/bar * bar = GHz

% CYH *** The following line is taken out, but **********************
% C       it is included in the ncpp term 
%c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
%      SUM = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      SUM = 0.0;

% CYH **************************************************************

      NLINES = length(F);
      for K = 1:NLINES
         DF = W300(K)*DEN;                  %Nico GHz/bar * bar = GHz
         FCEN = F(K);                       %Nico GHz
         Y = DEN * ( Y300(K) + V(K)*TH1 );  %Nico Y(T)

         STR = S300(K)*exp(-BE(K)*TH1);     %Nico [Hz*cm2]
         SF1 = (DF + (FREQ-FCEN)*Y)/((FREQ-FCEN)^2 + DF*DF); %Nico [1/GHz]
         SF2 = (DF - (FREQ+FCEN)*Y)/((FREQ+FCEN)^2 + DF*DF); %Nico [1/GHz]
         SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))^2;            %Nico [Hz*cm2/GHz] = [cm2]/1e9

      end

%      O2ABS = .5034E12*SUM*PRESDA*TH^3/3.14159;
%c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11      
%Nico The following computes n/pi*SUM*Th^2 (see notes)
%Nico n/pi = a/(pi*k*T_0) * Pda * Th
%Nico a/(pi*k*T_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont for units conversion (Pa->hPa, Hz->GHz, m->km)
%Nico Pa2hPa=1e-2; Hz2GHz=1e-9; m2cm=1e2; m2km=1e-3; Pa2hPa^-1 * Hz2GHz * m2cm^-2 * m2km^-1 = 1e-8
%Nico TH^3 = TH(from ideal gas law 2.13) * TH(from the MW approx of stimulated emission 2.16 vs. 2.14) * TH(from the partition sum 2.20)
      O2ABS = 1.6097E11*SUM*PRESDA*TH^3; 
      O2ABS = max(O2ABS,0.);

% CYH *** ********************************************************
% C   separate the equ. into line and continuum
% C   terms, and change the units from np/km to ppm

      npp = O2ABS;
%Nico intensities of the non-resonant transitions for O16-O16 and O16-O18, from JPL's line compilation
%c  1.571e-17 (o16-o16) + 1.3e-19 (o16-o18) = 1.584e-17
      ncpp = 1.584E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR));
%c   .20946e-4/(3.14159*1.38065e-19*300) = 1.6097e11      
%Nico a/(pi*k*T_0) = 0.20946/(3.14159*1.38065e-23*300) = 1.6097e19  - then it needs a factor 1e-8 to accont for units conversion (Pa->hPa, Hz->GHz)
%Nico Pa2hPa=1e-2; Hz2GHz=1e-9; m2cm=1e2; m2km=1e-3; Pa2hPa^-1 * Hz2GHz * m2cm^-2 * m2km^-1 = 1e-8
%Nico TH^3 = TH(from ideal gas law 2.13) * TH(from the MW approx of stimulated emission 2.16 vs. 2.14) * TH(from the partition sum 2.20)
      ncpp = 1.6097E11*ncpp*PRESDA*TH^3; %Nico: n/pi*SUM0
% C    add N2 term
      ncpp = ncpp + ABSN2_ros16(TEMP,PRES,FREQ);
% C     change the units from np/km to ppm
      npp = (npp /db2np)/factor;
      ncpp = (ncpp / db2np)/factor;

% CYH ************************************************************

return
      
end
      
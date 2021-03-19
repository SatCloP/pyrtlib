% C************************************
% C The interface of the following original function is changed to match
% C the interface of the ETL routine.    Yong Han, 1999.
% 
% C      FUNCTION ABH2O(T,P,RHO,F)
% C  Copyright (c) 2002 Massachusetts Institute of Technology
% C
% C  NAME- ABH2O_LIST    LANGUAGE- FORTRAN 77
% C
% C PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
% C 
%       IMPLICIT NONE
% C  CALLING SEQUENCE PARAMETERS-
% C    SPECIFICATIONS
%       REAL T,P,RHO,F,ABH2O
% C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
% C      T       KELVIN    I   TEMPERATURE
% C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
% C      RHO     G/M**3    I   WATER VAPOR DENSITY
% C      F       GHZ       I   FREQUENCY             
% C      ABH2O   NEPERS/KM O   POWER ABSORPTION COEFFICIENT
% C
% C   Multiply ABH2O by 4.343 to obtain dB/km.
% c
% c   Line parameters will be read from file h2o_list.asc; intensities should
% c   include the isotope abundance factors.
% C   This version uses a line-shape cutoff.
% C
% C   REVISION HISTORY-
% C    DATE- OCT.6, 1988 EQS AS PUBL.: P.W. Rosenkranz, CHAP. 2 in 
% C     ATMOSPHERIC REMOTE SENSING BY MICROWAVE RADIOMETRY 
% C     (M.A. Janssen, ed., 1993) (http://hdl.handle.net/1721.1/68611).
% C     OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
% C                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
% C     OCT. 24, 95  PWR -ADD 1 LINE.
% C     JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
% C                       REVISED CONTINUUM.
% C     Mar. 2, 2003   PWR - LINE SHIFT
% c     Nov. 3, 2012 intensities at base T=296K, get line param. from file.
% c     Aug. 6, 2015 read continuum param from the file also. 
% 
% CYH **** add the following lines ******************************
%
% REFERENCES FOR MEASUREMENTS (freq, Wair, Wself in GHZ/bar, Xair, Xself, D/Wair):
% (1) H. Liebe et al., JQSRT v.9, pp. 31-47 (1969)
% (2) V. Payne et al.,IEEE Trans. Geosci. Rem. Sens. v.46, pp.3601-3617 (2008)
% (3) G. Golubiatnikov, J. MOLEC. SPEC. vol. 230, pp.196-198 (2005)
% (4) M. Koshelev et al., J. Molec. Spec. v.241, pp.101-108 (2007)
% (5) J.-M. Colmont et al.,J. Molec. Spec. v.193, pp.233-243 (1999)
% (6) M. Tretyakov et al, JQSRT v.114 pp.109-121 (2013)
% (7) G. Golubiatnikov et al., JQSRT v.109, pp.1828-1833 (2008)
% (8) V. Podobedov et al., JQSRT v.87, pp. 377-385 (2004)
% (9) M. Koshelev, JQSRT v.112, pp.550-552 (2011)
% (10) D. Turner et al., IEEE Trans. Geosci. Rem. Sens. v.47 pp.3326-37 (2009),
%     re-adjusted for new line par. July 27, 2015
% Other parameters from HITRAN12.

%      subroutine H2O_xxx(pdrykpa,vx,ekpa,frq,npp,ncpp)
function [npp,ncpp] = h20_ros2016_xxx(pdrykpa,vx,ekpa,frq);

% Nico 2016/11/30 *********************************************************

% From h2o_lis.asc
%molecule freq,GHz   S(296K)    B      Wair  Xair  D/Wair  Wself  Xself      rotational Q.N.        Refs.
MTX = [...
  11   22.235080  0.1317E-13  2.144  2.665  0.76 -0.0088  13.60   1.00  6  1  6        5  2  3; ...%  1,2
  11  183.310087  0.2334E-11  0.668  2.936  0.77 -0.0240  14.76   0.85  3  1  3        2  2  0; ...%  2,3
  11  321.225630  0.7861E-13  6.179  2.426  0.67 -0.0590  10.65   0.54 10  2  9        9  3  6; ...%    4
  11  325.152888  0.2725E-11  1.541  2.847  0.64 -0.0045  13.95   0.74  5  1  5        4  2  2; ...%  4,5
  11  380.197353  0.2473E-10  1.048  2.831  0.54 -0.0278  14.40   0.89  4  1  4        3  2  1; ...%    6
  11  439.150807  0.2152E-11  3.595  2.024  0.63  0.0182   9.06   0.52  6  4  3        5  5  0; ...%    6
  11  443.018343  0.4494E-12  5.048  1.568  0.60  0.0000   7.96   0.50  7  5  2        6  6  1; ...%    6
  11  448.001085  0.2586E-10  1.405  2.587  0.66 -0.0464  13.01   0.67  4  2  3        3  3  0; ...%    6
  11  470.888999  0.8253E-12  3.597  2.153  0.66  0.0240   9.70   0.65  6  4  2        5  5  1; ...%    6
  11  474.689092  0.3274E-11  2.379  2.340  0.65 -0.0190  11.24   0.64  5  3  3        4  4  0; ...%    6
  11  488.490108  0.6721E-12  2.852  2.610  0.69  0.0690  13.58   0.72  6  2  4        7  1  7; ...%    6
  11  556.935985  0.1561E-08  0.159  3.115  0.69  0.0600  14.24   1.00  1  1  0        1  0  1; ...%    7
  11  620.700807  0.1704E-10  2.391  2.468  0.75  0.0000  11.94   0.68  5  3  2        4  4  1; ...%    8
  11  752.033113  0.1029E-08  0.396  3.114  0.68  0.0520  13.58   0.84  2  1  1        2  0  2; ...%    9
  11  916.171582   .4266E-10  1.441  2.698  0.72  -.0208  13.91   0.78  4  2  2        3  3  1];
% continuum terms
CTR = [300.          5.96E-10  3.        1.42E-8    7.5];                                          %    10

% From abh2o_list.f
% Read line parameters; units: GHz, Hz*cm^2, MHz/mb
REFTLINE = 296.0; % reference T for lines
FL = MTX(:,2);     % LINE FREQUENCIES GHz
S1 = MTX(:,3);     % LINE INTENSITIES AT 296K(?):
B2 = MTX(:,4);     % T COEFF. OF INTENSITIES:
W3 = MTX(:,5)/1e3; % AIR-BROADENED WIDTH PARAMETERS AT ???K (WAIR -> W3)
X  = MTX(:,6);     % T-EXPONENT OF AIR-BROADENING:
SR = MTX(:,7);     % RATIO OF SHIFT TO WIDTH
WS = MTX(:,8)/1e3; % SELF-BROADENED WIDTH PARAMETERS AT ???K (WSELF -> WS) 
XS = MTX(:,9);     % T-EXPONENT OF SELF-BROADENING:
% Read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON = CTR(1);
CF      = CTR(2);
XCF     = CTR(3);
CS      = CTR(4);
XCS     = CTR(5);

% Nico 2016/11/30 *********************************************************

% Just a trial to be removed! 
% if mimic_ros98 = 1; it changes parameters to first line (22.2 GHz) to the 98 values 
% if mimic_ros98 = 0; it is ininfluent 
mimic_ros98 = 0;
if mimic_ros98
REFTLINE = 300.0;   % accounts for -0.89 K
S1(1) = .1310E-13;  % accounts for -0.35 K
B2(1) = 2.144;      % no change
W3(1) = .00281;     % accounts for -4.00 K % it sccounts for the most (the value 0.002665 was suggested by Payne et al., 2008)
X(1)  = .69;        % accounts for +1.43 K
WS(1) = .01349;     % accounts for +0.12 K
XS(1) = .61;        % accounts for +0.39 K
CF    = 5.43E-10;   % accounts for +0.02 K
CS    = 1.8E-8;     % accounts for +0.22 K
end
%

%CYH ***********************************************
      db2np = log(10.) * 0.1;
      rvap = 0.01 * 8.314510 / 18.01528;
      factor = .182 * frq;
      T = 300./vx;
      P = (pdrykpa+ekpa)*10.;
      RHO = ekpa*10./(rvap*T);
      F = frq;
%CYH ***********************************************

      if RHO <= 0.
        ABH2O = 0.;
        npp = 0;
        ncpp = 0;
        return
      end
      PVAP = RHO * T / 217;
      PDA = P - PVAP;
      DEN = 3.335E16 * RHO; %! const includes isotopic abundance
%     CONTINUUM TERMS
      TI = REFTCON/T;
%     XCF and XCS include 3 for conv. to density & stimulated emission
      CON = ( CF * PDA * TI^XCF + CS * PVAP * TI^XCS ) * PVAP * F*F;

% Nico 2016/11/30 *********************************************************
%C      ADD RESONANCES
      NLINES = length(FL);
      TI = REFTLINE/T;
      TI2 = TI^2.5;
      SUM = 0.;
      for I = 1:NLINES
          WIDTHF = W3(I) * PDA * TI^X(I);
          WIDTHS = WS(I) * PVAP * TI^XS(I);
          WIDTH = WIDTHF + WIDTHS;
          SHIFT = SR(I)*WIDTHF;  % unknown temperature dependence for shift
          WSQ = WIDTH^2;
          %c  line intensities include isotopic abundance
          S = S1(I) * TI2 * exp(B2(I)*(1.-TI)); 
          DF(1) = F - FL(I) - SHIFT;
          DF(2) = F + FL(I) + SHIFT;
%C  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
          BASE = WIDTH/(562500. + WSQ);
%C  DO FOR POSITIVE AND NEGATIVE RESONANCES
          RES = 0.;
          for J = 1:2
              if abs(DF(J)) < 750. 
                  RES = RES + WIDTH/(DF(J)^2+WSQ) - BASE;
              end
          end
          SUM = SUM + S * RES * (F./FL(I))^2;
% Nico 2016/11/30 *********************************************************


%CYH **************************************************************
%C     separate the following original equ. into line and continuum
%C        terms, and change the units from np/km to ppm
%C      ABH2O = .3183E-4*DEN*SUM + CON
      npp = (.3183E-4 * DEN * SUM / db2np) / factor;
      ncpp = (CON / db2np) / factor;
%CYH *************************************************************
      end
      
return
end

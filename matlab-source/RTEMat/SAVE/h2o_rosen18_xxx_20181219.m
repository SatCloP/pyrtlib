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
% C This version should not be used with a line list older than June 2018,
% c nor the new list with an older version of this subroutine.
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
% c     June 19, 2018 changed file format, separate shift for self & foreign gas
% 
% CYH **** add the following lines ******************************
%
% REFERENCES FOR MEASUREMENTS (freq, Wair, Wself, Dair, Dself in GHZ/bar, Xair, Xself, XDair, XDself)
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
% (11) D. Turner et al., IEEE Trans. Geosci. Rem. Sens. v.47 pp.3326-37 (2009),
%    continuum re-adjusted for new line par. June 18, 2018
% Other parameters from HITRAN16. (updated June 14, 2018)
%
% Nico 2018/06/21 *********************************************************
% References:
%
% Rosenkranz, personal communication, 2018/06/20
% Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017

%      subroutine H2O_xxx(pdrykpa,vx,ekpa,frq,npp,ncpp)
function [npp,ncpp,SUM] = h20_ros2018_xxx(pdrykpa,vx,ekpa,frq);

% Nico 2018/02/18 *********************************************************

blk = -9999; % blank

% From h2o_lis.asc (2018)
% molecule freq,GHz  S(296K)    B     Wair  Xair  Wself Xself   Dair  XDair  Dself  XDself       Refs.
%          FL(i)     S1(i)    B2(i)   Wair  X(i)  Wself Xs(i)   Sair  Xh(I)  Sself  Xhs(I) - variable names in the code
MTX = [...
  11   22.235080  0.1335E-13  2.172  2.699  0.76  13.29  1.20  -0.033   2.6   0.814   blk ; ...% 1,2,10
  11  183.310087  0.2319E-11  0.677  2.945  0.77  14.78  0.78  -0.072   blk   0.173   blk ; ...% 2,3,10
  11  321.225630  0.7657E-13  6.262  2.426  0.73  10.65  0.54  -0.143   blk   0.278   blk ; ...% 4
  11  325.152888  0.2721E-11  1.561  2.847  0.64  13.95  0.74  -0.013   blk   1.325   blk ; ...% 4,5
  11  380.197353  0.2477E-10  1.062  2.868  0.54  14.40  0.89  -0.074   blk   0.240   blk ; ...% 6
  11  439.150807  0.2137E-11  3.643  2.055  0.69   9.06  0.52   0.051   blk   0.165   blk ; ...% 6
  11  443.018343  0.4440E-12  5.116  1.819  0.70   7.96  0.50   0.140   blk  -0.229   blk ; ...% 6
  11  448.001085  0.2588E-10  1.424  2.612  0.70  13.01  0.67  -0.116   blk  -0.615   blk ; ...% 6
  11  470.888999  0.8196E-12  3.645  2.169  0.73   9.70  0.65   0.061   blk  -0.465   blk ; ...% 6
  11  474.689092  0.3268E-11  2.411  2.366  0.71  11.24  0.64  -0.027   blk  -0.720   blk ; ...% 6
  11  488.490108  0.6628E-12  2.890  2.616  0.75  13.58  0.72  -0.065   blk  -0.360   blk ; ...% 6
  11  556.935985  0.1570E-08  0.161  3.115  0.75  14.24  1.00   0.187   blk  -1.693   blk ; ...% 7
  11  620.700807  0.1700E-10  2.423  2.468  0.79  11.94  0.75     0.0   blk   0.687  0.92 ; ...% 8 blanks in Dair are set to 0 (confirmed by Phil)
  11  752.033113  0.1035E-08  0.402  3.114  0.77  13.58  0.84   0.162   blk  -0.878   blk ; ...% 9
  11  916.171582  0.4275E-10  1.461  2.695  0.79  13.91  0.78     0.0   blk     0.0   blk];    %   blanks in Dair/Dself are set to 0 (confirmed by Phil)
% continuum terms
CTR = [300.          5.95E-10  3.        1.42E-8    7.5];                                    % 11

% From abh2o_list.f
% Read line parameters; units: GHz, Hz*cm^2, MHz/mb
REFTLINE = 296.0; % reference T for lines
FL = MTX(:,2);     % LINE FREQUENCIES [GHz]
S1 = MTX(:,3);     % LINE INTENSITIES AT 296K [Hz*cm^2]
B2 = MTX(:,4);     % T COEFF. OF INTENSITIES [unitless]
W3 = MTX(:,5)/1e3; % AIR-BROADENED WIDTH PARAMETERS AT REFTLINE [MHz/mb -> GHz/mb]
X  = MTX(:,6);     % T-EXPONENT OF AIR-BROADENING [unitless]
WS = MTX(:,7)/1e3; % SELF-BROADENED WIDTH PARAMETERS AT REFTLINE [MHz/mb -> GHz/mb]
XS = MTX(:,8);     % T-EXPONENT OF SELF-BROADENING [unitless]
SH = MTX(:,9)/1e3; % AIR-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
XH = MTX(:,10);    % T-EXPONENT OF AIR-SHIFTING [unitless]
SHS = MTX(:,11)/1e3; % SELF-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
XHS = MTX(:,12);   % T-EXPONENT OF SELF-SHIFTING [unitless]

% Replace non-existing shifting parameters with broadening parameters
indx = find(XH==blk); XH(indx) = X(indx);
indx = find(XHS==blk); XHS(indx) = XS(indx);

% Read continuum parameters; units: Kelvin, 1/(km*mb^2*GHz^2)
REFTCON = CTR(1);
CF      = CTR(2);
XCF     = CTR(3);
CS      = CTR(4);
XCS     = CTR(5);

% Nico 2018/02/18 *********************************************************

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
      PVAP = RHO * T / 216.68; % 1/rvap=216.67
      PDA = P - PVAP;
      DEN = 3.344E16 * RHO;  
      
%     CONTINUUM TERMS
      TI = REFTCON/T;
%     XCF and XCS include 3 for conv. to density & stimulated emission
      CON = ( CF * PDA * TI^XCF + CS * PVAP * TI^XCS ) * PVAP * F*F;

% Nico 2018/06/22 *********************************************************
%C      ADD RESONANCES
      NLINES = length(FL);
      TI = REFTLINE/T;
      TI2 = TI^2.5;
      SUM = 0.;
      for I = 1:NLINES
          WIDTHF = W3(I) * PDA * TI^X(I);
          WIDTHS = WS(I) * PVAP * TI^XS(I);
          WIDTH = WIDTHF + WIDTHS;
          SHIFTF = SH(I) * PDA * TI^XH(I);
          SHIFTS = SHS(I) * PVAP * TI^XHS(I);
          SHIFT = SHIFTF + SHIFTS;
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
          
% Nico 2018/06/22 *********************************************************


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

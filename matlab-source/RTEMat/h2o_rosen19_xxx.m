% Nico: from /Users/Nico/SFTWR/MPM/PWR2019/SpeedDependence_183GHz/new_abh2o/abh2o_list.f
% Nico: which, except for a comment, is the same as /Users/Nico/SFTWR/MPM/PWR2020/absorption/abh2o_list.f
%
% C************************************
% C The interface of the following original function is changed to match
% C the interface of the ETL routine.    Yong Han, 1999.
% 
%       FUNCTION ABH2O(T,P,RHO,F)
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
% C     Feb. 19, 2019 added Aair and Aself to file format
% 
% CYH **** add the following lines ******************************
%
% Nico 2019/03/19 *********************************************************
% References:
%
% Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017
% Rosenkranz, personal communication, 2019/03/18
% abh2o_list.f

%function [npp,ncpp,SUM] = h2o_rosen18_22sd(pdrykpa,vx,ekpa,frq);
function [npp,ncpp,SUM] = h2o_rosen19_xxx(pdrykpa,vx,ekpa,frq);

% Nico: the best-fit Voigt are given in Koshelev et al. 2018, Table 2 (RAD,
% MHz/Torr). These correspond to W3(1) and WS(1) in h2o_list_r19 (MHz/mb)

% Read the list of parameters
h2o_list_r19

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


% Nico 2019/03/18 *********************************************************
%C      ADD RESONANCES
      NLINES = length(FL);
      TI = REFTLINE/T;
      TILN = log(TI);
      TI2 = exp(2.5*TILN); % = TI^2.5;
      SUM = 0.;
      for I = 1:NLINES
          WIDTHF = WA(I) * PDA * TI^X(I);
          WIDTHS = WS(I) * PVAP * TI^XS(I);
          WIDTH = WIDTHF + WIDTHS;
          SHIFTF = SH(I) * PDA * (1. - Aair(I) * TILN ) * TI^XH(I);
          SHIFTS = SHS(I) * PVAP * (1. - Aself(I) * TILN ) * TI^XHS(I);        
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
                 RES = RES + WIDTH / ( DF(J)^2 + WSQ ) - BASE;
              end
          end
          SUM = SUM + S * RES * (F./FL(I))^2;
          
% Nico 2019/03/19 *********************************************************

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

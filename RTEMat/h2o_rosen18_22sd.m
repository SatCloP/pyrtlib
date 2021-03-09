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
% c     Dec. 11, 2018 speed-dependent parameters only for 22 Ghz line
% 
% CYH **** add the following lines ******************************
%
% Nico 2018/12/19 *********************************************************
% References:
%
% Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017
% Rosenkranz, personal communication, 2018/12/13
% abh2o_22sd.f

function [npp,ncpp,SUM] = h2o_rosen18_22sd(pdrykpa,vx,ekpa,frq);

%c  correction factors (rel. to best-fit Voigt) for speed-dependence at 22 GHz:
R = [ 1.0153 .1612 1.0258 .1440]; % r0air, r2air, r0self, r2self
% Nico: the best-fit Voigt are given in Koshelev et al. 2018, Table 2 (RAD,
% MHz/Torr). These correspond to W3(1) and WS(1) in h2o_list_r18 (MHz/mb)

% Read the list of parameters
h2o_list_r18

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

% Nico 2018/12/19 *********************************************************
%C      ADD RESONANCES
      NLINES = length(FL);
      TI = REFTLINE/T;
      TI2 = TI^2.5;
      SUM = 0.;
      for I = 1:NLINES
          WIDTHF = W3(I) * PDA * TI^X(I);
          WIDTHS = WS(I) * PVAP * TI^XS(I);
          if I == 1 % 22 GHz
             WIDTH0 = WIDTHF*R(1) + WIDTHS*R(3);
             WIDTH2 = WIDTHF*R(2) + WIDTHS*R(4);
          else
             WIDTH0 = WIDTHF + WIDTHS;
             WIDTH2 = 0.;
          end
          SHIFTF = SH(I) * PDA * TI^XH(I);
          SHIFTS = SHS(I) * PVAP * TI^XHS(I);
          SHIFT = SHIFTF + SHIFTS;
          %c SHIFT2 is too small to measure
          % Nico: thus using the best-fit Voigt (SHIFT instead of SHIFT0 and SHIFT2)
          WSQ = WIDTH0^2;
          
          %c  line intensities include isotopic abundance
          S = S1(I) * TI2 * exp(B2(I)*(1.-TI)); 
          DF(1) = F - FL(I) - SHIFT;
          DF(2) = F + FL(I) + SHIFT;
          %C  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
          BASE = WIDTH0/(562500. + WSQ);
          %C  DO FOR POSITIVE AND NEGATIVE RESONANCES
          RES = 0.;
          for J = 1:2
             %IF(I.EQ.1 .AND. J.EQ.1 .AND. ABS(DF(J)).LT.10.*WIDTH0) THEN
             if I==1 & J==1 & abs(DF(J)) < 10*WIDTH0 
                %C speed-dependent resonant shape factor
                % double complex dcerror,Xc,Xrt,pxw,A
                Xc = complex((WIDTH0-1.5*WIDTH2),DF(J))/WIDTH2;
                Xrt = sqrt(Xc);                
                %pxw = 1.77245385090551603 * Xrt * dcerror(-imag(Xrt),dble(Xrt));
                pxw = 1.77245385090551603 * Xrt * dcerror(-imag(Xrt),real(Xrt)); % ATT!! I'm not totally sure dble correspond to real
                A = 2. * (1. - pxw) / WIDTH2;                
                RES = RES + real(A) - BASE;                 
             elseif abs(DF(J)) < 750. 
                RES = RES + WIDTH0/(DF(J)^2+WSQ) - BASE;
             end
          end
          SUM = SUM + S * RES * (F./FL(I))^2;
          
% Nico 2018/12/19 *********************************************************

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

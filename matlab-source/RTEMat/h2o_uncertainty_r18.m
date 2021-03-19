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
%      continuum re-adjusted for new line par. July 16, 2018.
% Other parameters from HITRAN2016. (updated July 16, 2018)
%
% Nico 2018/06/21 *********************************************************
% References:
%
% Rosenkranz, personal communication, 2018/06/20
% Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017

%      subroutine H2O_xxx(pdrykpa,vx,ekpa,frq,npp,ncpp)
function [npp,ncpp,SUM] = h2o_uncertainty_r18(pdrykpa,vx,ekpa,frq,AMU);

% Read the list of parameters
h2o_list_r18

% Perturb parameters ******************************************************
CF = AMU.con_Cf.value*AMU.con_Cf_factr.value; % foreign continuum; only the factor is perturbed
CS = AMU.con_Cs.value*AMU.con_Cs_factr.value; % self continuum; only the factor is perturbed
XCF = AMU.con_Xf.value; % temp dependence fore cont; 
XCS = AMU.con_Xs.value; % temp dependence self cont; 
S1(1:2) = AMU.S.value(1:2); % line intensities
B2 = AMU.B2.value; % Temp. coeff. of intensity
W3(1:2) = AMU.gamma_a.value(1:2) / 1e3; % air-pressure broadening
WS(1:2) = AMU.gamma_w.value(1:2) / 1e3; % self broadening
X(1:2) = AMU.n_a.value(1:2); % temp exponent of air broadening
XS(1:2) = AMU.n_w.value(1:2); % temp exponent of self broadening
%SR(1:2) = AMU.SR.value(1:2); % shift to width ratio
FL(1:2) = AMU.FL.value(1:2); % line frequency
n_S = AMU.wv_nS.value; % intensity temperature-dependence exponent
SH(1:2) = AMU.delta_a.value(1:2) / 1e3; % AIR-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
SHS(1:2) = AMU.delta_w.value(1:2) / 1e3; % SELF-BROADENED SHIFT PARAMETERS [MHz/mb -> GHz/mb]
%**************************************************************************
% Nico 2018/08/08 *********************************************************

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
      %TI2 = TI^2.5;
      TI2 = TI^n_S;
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

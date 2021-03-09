% C     COMPUTES POWER ABSORPTION COEFFICIENT IN NEPERS/KM 
% C     BY SUSPENDED CLOUD LIQUID WATER DROPLETS. MULTIPLY ABLIQ BY
% C     4.343 TO CONVERT TO DB/KM.
% C
% c     ARGUMENTS (INPUT):
% C     WATER IN G/M**3
% C     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
% C     TEMP IN KELVIN
% C
% C     REVISION HISTORY:
% C        PWR 6/5/15   using dilec12 for complex dielectric constant
% C
%
% Nico - 2018/02/19
%
% References:
%
% Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017

%      FUNCTION ABLIQ(WATER,FREQ,TEMP)
function ABLIQ = ABLIQ(WATER,FREQ,TEMP);

      if WATER <= 0
         ABLIQ = 0;
         return
      end
      
      EPS = dilec12(FREQ,TEMP);
      
%       THETA1 = 1. - 300. / TEMP;
%       EPS0 = 77.66 - 103.3 * THETA1;
%       EPS1 = .0671 * EPS0;
%       EPS2 = 3.52;                  % from MPM93
%       FP = 20.1 * exp(7.88*THETA1); % from eq. 2b
%       FS = 39.8 * FP;
%       EPS = (EPS0 - EPS1) / complex(1.,FREQ/FP) + (EPS1-EPS2) / complex(1.,FREQ/FS) + EPS2;

      RE = (EPS - 1.) / (EPS + 2.);
      ABLIQ = -0.06286 * imag(RE) * FREQ * WATER;
      
      
return
      
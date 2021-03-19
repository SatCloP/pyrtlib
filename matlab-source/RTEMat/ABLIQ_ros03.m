%     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
%     ARGUMENTS (INPUT):
%     WATER IN G/M**3
%     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
%     TEMP IN KELVIN
%
%     REFERENCES:
%     LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
%      (1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.
%
%     REVISION HISTORY:
%        PWR 8/3/92   original version
%        PWR 12/14/98 temp. dependence of EPS2 eliminated to agree 
%                     with MPM93 
%        pwr 2/27/02  use exponential dep. on T, eq. 2b instead of eq. 4a 
%

%      FUNCTION ABLIQ(WATER,FREQ,TEMP)
function ABLIQ = ABLIQ(WATER,FREQ,TEMP);

      if WATER <= 0
         ABLIQ = 0;
         return
      end
      
      THETA1 = 1. - 300. / TEMP;
      EPS0 = 77.66 - 103.3 * THETA1;
      EPS1 = .0671 * EPS0;
      EPS2 = 3.52;                  % from MPM93
      FP = 20.1 * exp(7.88*THETA1); % from eq. 2b
      FS = 39.8 * FP;
      EPS = (EPS0 - EPS1) / complex(1.,FREQ/FP) + (EPS1-EPS2) / complex(1.,FREQ/FS) + EPS2;
      RE = (EPS - 1.) / (EPS + 2.);
      ABLIQ = -0.06286 * imag(RE) * FREQ * WATER;
      
return
      
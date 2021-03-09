%       FUNCTION ABLIQ(WATER,FREQ,TEMP)
% C     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
% c     ARGUMENTS (INPUT):
% C     WATER IN G/M**3
% C     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
% C     TEMP IN KELVIN
% C
% C     REFERENCES:
% C     LIEBE, HUFFORD AND MANABE, INT. J. IR & MM WAVES V.12, pp.659-675
% C      (1991);  Liebe et al, AGARD Conf. Proc. 542, May 1993.
% c
% C     REVISION HISTORY:
% C        PWR 8/3/92   original version
% c        PWR 12/14/98 temp. dependence of EPS2 eliminated to agree 
% c                     with MPM93 
% C

function ABLIQ = ABLIQ(WATER,FREQ,TEMP);

      if WATER <= 0
         ABLIQ = 0;
         return
      end
 
      THETA1 = 1. - 300. / TEMP;
      EPS0 = 77.66 - 103.3 * THETA1;
      EPS1 = .0671 * EPS0;
      EPS2 = 3.52;
      FP = (316. * THETA1 + 146.4) * THETA1 + 20.20;
      FS = 39.8 * FP;
      EPS = (EPS0 - EPS1) / complex(1.,FREQ/FP) + (EPS1-EPS2) / complex(1.,FREQ/FS) + EPS2;
      RE = (EPS - 1.) / (EPS + 2.);
      ABLIQ = -0.06286 * imag(RE) * FREQ * WATER;
      
return
      
%       FUNCTION ABSN2(T,P,F)
% C  Copyright (c) 2002 Massachusetts Institute of Technology
% C     ABSN2 = COLLISION-INDUCED POWER ABSORPTION COEFFICIENT 
% C     (NEPER/KM) IN AIR
% C     T = TEMPERATURE (K)
% C     P = DRY AIR PRESSURE (MB)
% C     F = FREQUENCY (GHZ)(valid 0-2000 GHz)
% C
% c     5/22/02, 4/14/05 P.Rosenkranz
% c
% c     References:
% c     See eq. 2.6 in Thermal Microwave Radiation - Applications 
% c      for Remote Sensing (C. Maetzler, ed.) London, IET, 2006.
% C     Based on model by:
% C      Borysow, A, and L. Frommhold, 
% C      Astrophysical Journal, v.311, pp.1043-1057 (1986)
% C     with modification of 1.34 to account for O2-O2 and O2-N2
% c     collisions, as calculated by
% c      J.Boissoles, C.Boulet, R.H.Tipping, A.Brown, and Q.Ma, 
% C      J. Quant. Spectros. Radiat. Trans. v.82, 505-516 (2003).
% c

function ABSN2 = ABSN2(T,P,F);

      TH = 300./T;
      FDEPEN = .5 + .5/(1.+(F/450.).^2);
      BF = 6.5E-14*FDEPEN*P*P*F*F*TH^3.6;
      ABSN2 = 1.34*BF;

%       RETURN
%       END

return

%       FUNCTION ABSN2(T,P,F)
% C  Copyright (c) 2002 Massachusetts Institute of Technology
% C     ABSN2 = COLLISION-INDUCED POWER ABSORPTION COEFFICIENT 
% C     (NEPER/KM) IN AIR ("dry continuum", mostly due to N2-N2, 
% C     but also contributions from O2-N2 and O2-O2)
% C     T = TEMPERATURE (K)
% C     P = DRY AIR PRESSURE (MB)
% C     F = FREQUENCY (GHZ)(valid 0-2000 GHz)
% C       Multiply ABSN2 by 4.343 to obtain dB/km.
% C
% c     5/22/02, 4/14/05, 6/23/18 P.Rosenkranz
% c
% c     References:
% C     Frequency dependence based on model by A. Borysow and L. Frommhold,
% C      Astrophysical Journal, v.311, pp.1043-1057 (1986).
% c     See Eq. 2.6 in Thermal Microwave Radiation - Applications 
% c      for Remote Sensing (C. Maetzler, ed.) London, IET, 2006.
% c     Amplitude increased by 14% based on analysis by
% c      M. Tretyakov and A. Zibarova, JQSRT v.216, pp. 70-75 (2018)  
% c

function ABSN2 = ABSN2(T,P,F);

      TH = 300./T;
      FDEPEN = .5 + .5/(1.+(F/450.)^2);
      ABSN2 = 9.95E-14*FDEPEN*P*P*F*F*TH^3.22;
%      RETURN
%      END

return

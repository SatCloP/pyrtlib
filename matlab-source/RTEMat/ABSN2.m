%       FUNCTION ABSN2(T,P,F)
% C     ABSN2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
% C             (NEPER/KM)
% C     T = TEMPERATURE (K)
% C     P = PRESSURE (MB)
% C     F = FREQUENCY (GHZ)
% C
function ABSN2 = ABSN2(T,P,F);

      TH = 300./T;
      ABSN2 = 6.4E-14*P*P*F*F*TH^3.55;

%       RETURN
%       END

return

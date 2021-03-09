% C  SIXTH-ORDER APPROX TO THE COMPLEX ERROR FUNCTION OF z=X+iY.
% C  CERROR = exp(-z^2)erfc(-iz)
% C  This version is double precision and valid in all quadrants.
% C  REFERENCE FOR EQUATIONS:
% C      HUI, ARMSTRONG AND WRAY, JQSRT V.19, P.509-516 (1978).
% C
% C     P. ROSENKRANZ  12/11/2018
%
% 2018/12/19 - Nico: first created from dcerror.f

function DCERROR = dcerror(X,Y)

%      IMPLICIT NONE
%      DOUBLE PRECISION X,Y,a(0:6),b(0:6)
%      DOUBLE COMPLEX ASUM,BSUM,ZH,w

a = [ 122.607931777104326, 214.382388694706425,...
      181.928533092181549,  93.155580458138441,...
       30.180142196210589,   5.912626209773153,...
        0.564189583562615];
    
b = [ 122.607931773875350, 352.730625110963558,...
      457.334478783897737, 348.703917719495792,...
      170.354001821091472,  53.992906912940207,...
       10.479857114260399];

% c  compute w in quadrants 1 or 2
% c  from eqs.(13), w(z) = [w(-z*)]*
% c  expansion in terms of ZH results in conjugation of w when X changes sign.
ZH = complex(abs(Y),-X);
ASUM = ((((( a(7)*ZH + a(6))*ZH + a(5))*ZH + a(4))*ZH + a(3))*ZH + a(2))*ZH + a(1);
BSUM = (((((( ZH+b(7)) * ZH+b(6)) * ZH+b(5)) * ZH+b(4)) * ZH+b(3))*ZH+ b(2))*ZH + b(1);
w = ASUM/BSUM;

if Y >= 0
   DCERROR = w; % for quadrants 1&2
else
   %c from eqs.(13), w(z) = 2exp(-z^2)-[w(z*)]*
   DCERROR = 2. * exp(-complex(X,Y)^2) - conj(w); % for quadrants 3&4
end

return
end

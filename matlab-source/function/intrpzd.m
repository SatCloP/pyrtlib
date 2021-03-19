% THIS FUNCTION COMPUTES THE INTEGRAL OF A DISCRETE FUNCTION.
% IT APPLIES SIMPLE QUADRATURE (TRAPEZOIDAL ALGORITHM).
% IT WANTS IN INPUT THE DIPENDENT AND INDIPENDENT VECTOR.
% ES: [INT]=intrpzd(x,y)
%
% N.B.: Already exist a Matlab function "trapz"
%       I need to check compatibility between these two.
%
% Nico, 2000
%
% 2013/05/21 - Corrected a bug (if x,y are empty it gives NaN) 

function [INT]=intrpzd(x,F)

if isempty(x) | isempty(F)
   INT = NaN;
   return
end

INT=0;

for i=1:length(x)-1
         
   Dx=x(i+1)-x(i);
   MF=(F(i+1)+F(i))/2;
   INT=INT+MF*Dx;
   
end

return
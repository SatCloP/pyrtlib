% findjumps introduces NaNs in a time series when there  
% are time jumps larger than a threshold given in input.
%
% [Xn,Yn,NofNaNs]=findjumps(X,Y,THR)

function [Xn,Yn,NofNaNs]=findjumps(X,Y,THR)

if length(X) ~= length(Y)
   fprintf('findjumps: Different size for X and Y.\n')
   return
end

DX = diff(X);
indx = find(DX > THR);
flag = zeros(length(X),1);
flag(indx+1) = 1;

Xn = zeros(length(X)+length(indx),1);
Yn = zeros(length(X)+length(indx),1);

in=0;
for i=1:length(X)
   
   if flag(i) 
      Xn(i+in)=NaN;
      Yn(i+in)=NaN;
      in=in+1;
   end
   
   Xn(i+in) = X(i);
   Yn(i+in) = Y(i);
   
end

NofNaNs = in;

return

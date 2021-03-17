function psat = eswat(T);

% function psat = eswat(T);
%
% determine saturation vapor pressure of H2O w/r/t liquid water
% T in tmperature in degrees K
% This is taken from Barb W. from a pamphlet of notes originally 
% from Steven Ackerman
%
% psat is the saturation pressure in mbar

A = [...
  6984.505294 
 -188.9039310 
  2.133357675
 -1.288580973E-2
  4.393587233E-5
 -8.023923082E-8
  6.136820929E-11];

psat = A(1)+T.*(A(2)+T.*(A(3)+T.*(A(4)+T.*(A(5)+T.*(A(6)+A(7).*T)))));

ind = find(psat < 0 | T < 200);
psat(ind) = psat(ind).*0;


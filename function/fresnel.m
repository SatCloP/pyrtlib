% Function fresnel
%
% This function calculates the Horizontal and the Vertical reflectivity
% for sea water (?every media?) knowing the complex dielectric constant epsilon
% and the observation angle with the surface normal (degrees).
%
% Es:
%    [rH,rV]=fresnel(eps,tht)
%
% INPUTS:
%        eps: complex dielectric constant (as calculated by dielectriconstant.m)
%        tht: observation angle with the surface normal (degrees)
%
% NOTES:
%        The function [eps,rH,rV]=epsilon(sol,temp,wave,tht) 
%        in C:/Nico/5mm/mat5mm is now obsolete
%
% Nico, 8/2000

function [rH,rV]=fresnel(eps,tht)

  d2r=pi/180;
  tht=tht*d2r;
  
  C = (eps-sin(tht)^2)^(1/2);
  CT= cos(tht);
  CRH = (CT-C)/(CT+C);
  CRV = (eps*CT-C)/(eps*CT+C);
  
  %CTH=2*CT/(CT+C);
  %CTV=2*eps*CT/(eps*CT+C);
  
  rH=real(CRH*conj(CRH));  % It's the same of ( abs(CRH) )^2
  rV=real(CRV*conj(CRV));  % It's the same of ( abs(CRV) )^2

return

function radiance = bt2rad(freq,bt);

%
% function radiance = bt2rad(freq,bt);
%
% compute radiance given wavenumbers and brightness temperature.
%
% Inputs:
%	freq	  wavenumbers (Nx1) in 1/cm 
%	bt        brightnes temperature (Nx1) in Kelvin
% Outputs:
%	radiance  Planck radiances (Nx1) in mW/(m^2 sr. 1/cm)
%
% see also RADTOT.M, RAD2BT.M, TTORAD.M
% 
% DCT 11/11-99
%

% fundamental constants:
%  (Cohen, E.R. and B.N. Taylor, The 1986 CODATA recommended values
%  of the fundamental physical constants, Journal of Research of
%  the National Bureau of Standard, 92(2), March-April 1987.)

h = 6.6260755E-34;  % Planck constant in Js
c = 2.99792458E8;   % photon speed in m/s
k = 1.380658E-23;   % Boltzmann constant in J/K

c1 = 2*h*c*c*1e8;
c2 = h*c/k*1e2;

radiance = 1e3 * c1.*(freq.*freq.*freq)./(exp((c2.*freq)./bt)-1);

ind = find(bt <= 0 | ~isreal(bt));
if length(ind) ~= 0
   disp('WARNING: negative and/or imaginary inputs to BT2RAD')
end

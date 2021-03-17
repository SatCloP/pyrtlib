function [e1,e2] = rh2e(p,t,rh,Tconvert);

% function [e1,e2] = rh2e(p,t,rh,Tconvert);
%
% determine H2O partial pressure (e,mbar) given
% pressure (p,mbar), temperature (t,K), and
% relative humidity (rh,%)
%
% Two H2O partial pressures are returned: e1 is with RH defined as 
% the ratio of water vapor partial pressure to saturation vapor
% pressure and e2 is with RH defined as the ratio of water vapor 
% mixing ratio to saturation mixing ratio.
%
% if input, Tconvert is used as the temperature point to switch
% from using saturation vapor pressure over water to over ice.
%
% DCT 3/6/00
%

% ratio of water mass to dry air mass
eps = 0.621970585;

% saturation pressure
if nargin == 3
  esat = satvap(t); % Goff Gratch formulation, over water
  wsat = satmix(p,t); % Goff Gratch formulation, over water
elseif nargin == 4
  esat = satvap(t,Tconvert); % Goff Gratch formulation, over water/ice
  wsat = satmix(p,t,Tconvert); % Goff Gratch formulation, over water/ice
end

% H2O partial pressure w/ RH defined as ratio of pressures
% This is the one used within RTE; see Vapor_xxx.m : e = rh .* es;
e1 = rh./100.*esat;

% H2O partial pressure w/ RH defined as ratio of mixing ratios
% rh = 100*w./wsat = 100*(e/(p-e))/(esat/(p-esat))
e2 = (rh/100.*esat./(p-esat).*p)./(rh/100.*esat./(p-esat) + 1);


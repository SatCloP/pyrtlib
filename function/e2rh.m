function [rh1,rh2] = e2rh(p,t,e,Tconvert);

% function [rh1,rh2] = e2rh(p,t,e,Tconvert);
%
% determine relative humidity (e,mbar) given
% pressure (p,mbar), temperature (t,K), and
% H2O partial pressure (e,mbar)
%
% Two RHs are returned: rh1 is with RH defined as the ratio
% of water vapor partial pressure to saturation vapor pressure and
% rh2 is with RH defined as the ratio of water vapor mixing ratio to 
% saturation mixing ratio.
%
% if input, Tconvert is used as the temperature point to switch
% from using saturation vapor pressure over water to over ice.
%
% DCT 3/6/00
%

% saturation pressure
if nargin == 3
  esat = satvap(t); % Goff Gratch formulation, over water
  wsat = satmix(p,t); % Goff Gratch formulation, over water
elseif nargin == 4
  esat = satvap(t,Tconvert); % Goff Gratch formulation, over water/ice
  wsat = satmix(p,t,Tconvert); % Goff Gratch formulation, over water/ice
end

% w/ RH defined as ratio of pressures
rh1 = 100.*e./esat;

% w/ RH defined as a ratio of mass mixing ratios
w = e2mr(p,e);
rh2 = 100.*w./wsat;

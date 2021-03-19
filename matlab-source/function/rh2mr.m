function [w1,w2] = rh2mr(p,t,rh,Tconvert);

% function [w1,w2] = rh2mr(p,t,rh,Tconvert);
%
% determine H2O mixing ratio (w, g/kg) given
% reference pressure (mbar), temperature (t,K), and
% relative humidity (rh,%)
%
% Two mixing ratios are returned: w1 is with RH defined as the ratio 
% of water vapor partial pressure to saturation vapor pressure and
% w2 is with RH defined as the ratio of water vapor mixing ratio to 
% saturation mixing ratio.
%
% if input, Tconvert is used as the temperature point to switch
% from using saturation vapor pressure over water to over ice.
%
% DCT 3/5/00
%

% saturation pressure
if nargin == 3
  esat = satvap(t); % Goff Gratch formulation, over water
  wsat = satmix(p,t); % Goff Gratch formulation, over water
elseif nargin == 4
  esat = satvap(t,Tconvert); % Goff Gratch formulation, over water/ice
  wsat = satmix(p,t,Tconvert); % Goff Gratch formulation, over water/ice
end

% H2O partial pressure
e = rh./100.*esat;

% H2O mixing ratio
w1 = e2mr(p,e);

% using WMO definition of relative humidity
w2 = rh./100.*wsat;


function [rh1,rh2] = mr2rh(p,t,w,Tconvert);

%
% function [rh1,rh2] = mr2rh(p,t,w,Tconvert);
%
% determine relative humidity (%) given
% reference pressure (mbar), temperature (t,K), and
% water vapor mass mixing ratio (w,g/kg)
%
% Two RHs are returned: rh1 is with RH defined as the ratio 
% of water vapor partial pressure to saturation vapor pressure and
% rh2 is with RH defined as the ratio of water vapor mixing ratio to 
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
e = mr2e(p,w);

% RH using ratios of gas pressures
rh1 = 100.*e./esat;

% RH using WMO definition of relative humidity
rh2 = 100*w./wsat;


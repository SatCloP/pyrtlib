
function [rh1,rh2] = dp2rh(p,t,dp,Tconvert)

%
% function [rh1,rh2] = dp2rh(p,t,dp,Tconvert)
%
% compute relative humidity (%) given total pressure p (mb), 
% air temperature t (K), and dew point temperature (K).
%
% if input, Tconvert is used as the AIR temperature to switch
% from using saturation vapor pressure over water to over ice.
%
% Two relative humidities are returned: rh1 is with RH defined as
% the ratio of water vapor partial pressure to saturation vapor
% pressure and rh2 is with RH defined as the ratio of water vapor 
% mixing ratio to saturation mixing ratio.
%
% dct 3/6/2000
%

% vapor pressures computed at dew point temperature over water and ice
ewater = eswat_goffgratch(dp); % Goff Gratch formulation, over water
eice = esice_goffgratch(dp); % Goff Gratch formulation, over ice

% saturation pressure and mixing ratios computed at air temperaure
if nargin == 3
  esat = satvap(t);
  wsat = satmix(p,t);
elseif nargin == 4
  esat = satvap(t,Tconvert);
  wsat = satmix(p,t,Tconvert);
end

e = ewater;  % initialize water vapor partial pressure 
             % as over water for all air temperatures

if nargin == 4
  ind = find(t <= Tconvert);
  e(ind) = eice(ind);   % over ice for air temperatures <= Tconvert
end

% water vapor mixing ratio
w = e2mr(p,e);

% corresponding RH (as a ratio of pressures)
rh1 = 100.*e./esat;

% WMO definition of RH gives a second mixing ratio and RH
rh2 = 100.*w./wsat;


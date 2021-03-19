function w = dp2mr(p,t,dp,Tconvert)

%
% function w = dp2mr(p,t,dp,Tconvert)
%
% compute water vapor mixing ratio (g/kg) given total 
% pressure p (mb), air temperature t (K), and dew point 
% temperature (K).
%
% if input, Tconvert is used as the AIR temperature to switch
% from using saturation vapor pressure over water to over ice.
%
% dct 3/5/2000
%

% vapor pressures computed at dew point temperature over water and ice
ewater = eswat_goffgratch(dp); % Goff Gratch formulation, over water
eice = esice_goffgratch(dp); % Goff Gratch formulation, over ice

% initialize water vapor partial pressure as over water
e = ewater;
% use over ice for air temperatures <= Tconvert
if nargin == 4
  ind = find(t <= Tconvert);
  e(ind) = eice(ind);
end

% water vapor mixing ratio
w = e2mr(p,e);

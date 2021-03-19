
function e = dp2e(t,dp,Tconvert)

%
% function e = dp2e(t,dp,Tconvert)
%
% compute H2O partial pressure (e,mabr) given total pressure p (mb), 
% air temperature t (K), and dew point temperature (K).
%
% if input, Tconvert is used as the AIR temperature to switch
% from using saturation vapor pressure over water to over ice.
%
% dct 3/6/2000
%

% vapor pressures computed at dew point temperature over water and ice
ewater = eswat_goffgratch(dp); % Goff Gratch formulation, over water
eice = esice_goffgratch(dp); % Goff Gratch formulation, over ice

e = ewater;  % initialize water vapor partial pressure 
             % as over water for all air temperatures

if nargin == 3
  ind = find(t <= Tconvert);
  e(ind) = eice(ind);   % over ice for air temperatures <= Tconvert
end

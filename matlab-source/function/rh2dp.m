
function [dp1,dp2] = rh2dp(p,t,rh,Tconvert)

%
% function [dp1,dp2] = rh2dp(p,t,rh,Tconvert)
%
% Compute dew point temperature given relative humidity (%).  
% Uses Goff-Gratch formulation.
%
% If input, the air temperature (t,K) is used to determine if
% saturation over water (for t > Tconvert) or saturation over ice
% (for t <= Tconvert) is used.
%
% Two dew points are returned: dp1 is with RH defined as the ratio 
% of water vapor partial pressure to saturation vapor pressure and
% dp2 is with RH defined as the ratio of water vapor mixing ratio to 
% saturation mixing ratio.
%
% Notes: results not valid for dew points >= 370 K and  <= 160 K.
%
% Reference: Goff-Gratch formulation from sixth revised 
%       edition of Smithsonian Meteorology Tables.
%
% DCT 3/6/01
%

% compute H2O partial pressures using traditional and WMO
% definitions of RH
if nargin < 4
  [e1,e2] = rh2e(p,t,rh);
elseif nargin == 4
  [e1,e2] = rh2e(p,t,rh,Tconvert);
end

% interpolate (saturation pressure vs T) to desired pressure
T = (100:0.2:400)';
esat_water = eswat_goffgratch(T);
dp1 = interp1(esat_water,T,e1);
dp2 = interp1(esat_water,T,e2);

% over ice for air temperature <= Tconvert
if nargin == 4
  esat_ice = esice_goffgratch(T);
  ind = find(t <= Tconvert);
  dp1(ind) = interp1(esat_ice,T,e1(ind));
  dp2(ind) = interp1(esat_ice,T,e2(ind));
end


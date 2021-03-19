
function dp = e2dp(e,t,Tconvert)

%
% function dp = e2dp(e,t,Tconvert)
%
% Compute dew point temperature given water vapor saturation 
% pressure (e, mbar).  Uses Goff-Gratch formulation.
%
% If input, the air temperature (t,K) is used to determine if
% saturation over water (for t > Tconvert) or saturation over ice
% (for t <= Tconvert) is used.
%
% Notes: results not valid for dew points >= 370 K and  <= 160 K.
%
% Reference: Goff-Gratch formulation from sixth revised 
%       edition of Smithsonian Meteorology Tables.
%
% DCT 3/6/01
%

T = (100:0.2:400)';
esat_water = eswat_goffgratch(T);
dp = interp1(esat_water,T,e);

if nargin == 3
  esat_ice = esice_goffgratch(T);
  ind = find(t <= Tconvert);
  dp(ind) = interp1(esat_ice,T,e(ind));
end


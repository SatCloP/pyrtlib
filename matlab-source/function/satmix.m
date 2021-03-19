
function wsat = satmix(p,T,Tconvert)

%
% function wsat = satmix(p,T,Tconvert)
%
% compute saturation mixing ratio [g/kg] given reference pressure, 
% p [mbar] and temperature, T [K].  If Tconvert input, the calculation uses 
% the saturation vapor pressure over ice (opposed to over water) 
% for temperatures less than Tconvert [K].
%
% DCT, updated 3/5/00
%

warning off

% saturation pressure
if nargin == 2
  esat = satvap(T);
elseif nargin == 3
  esat = satvap(T,Tconvert);
end

% saturation mixing ratio
wsat = e2mr(p,esat);

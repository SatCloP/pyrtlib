function esat = satvap(T,Tconvert)

% function esat = satvap(T,Tconvert)
%
% compute saturation vapor pressure [mbar] given temperature, T [K].
% If Tconvert is input, the calculation uses the saturation vapor 
% pressure over ice (opposed to over water) for temperatures less than 
% Tconvert [K].
%
% DCT, updated 3/5/00
%

% saturation pressure over water
esat = eswat_goffgratch(T); % Goff Gratch formulation, over water

% saturation pressure over ice if needed
if nargin == 2
  ind = find(T <= Tconvert);
  esat(ind) = esice_goffgratch(T(ind)); % Goff Gratch formulation, over ice
end


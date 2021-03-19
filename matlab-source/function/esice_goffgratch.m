function svp = esice_goffgratch(T)

%
% function svp = esice_goffgratch(T)
%
% Compute water vapor saturation pressure over ice
% using Goff-Gratch formulation.  Adopted from PvD's 
% svp_ice.pro.
%
% Inputs: 
%    T       temperature [Kelvin]
% 
% Output:
%    svp    saturation pressure [mbar]
%
% Notes: svp returned for all values of input T,
%    but results not valid for T >= 370 K and 
%    T <= 160 K.
%
% Reference: Goff-Gratch formulation from sixth revised 
%       edition of Smithsonian Meteorology Tables.
%
% DCT 8/22/00
%

ewi = 6.1071;
c1 =  9.09718;
c2 = 3.56654;
c3 = 0.876793;

ratio = 273.15 ./ T;

tmp = ( -c1 * ( ratio - 1.0 ) ) - ...
      (  c2 * log10( ratio ) ) + ...
      (  c3 * ( 1.0 - ( 1.0 ./ ratio ) ) ) + ...
      log10( ewi );

svp = 10.0.^tmp;


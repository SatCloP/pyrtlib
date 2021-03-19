function svp = eswat_goffgratch(T)

%
% function svp = eswat_goffgratch(T)
%
% Compute water vapor saturation pressure over water
% using Goff-Gratch formulation.  Adopted from PvD's 
% svp_water.pro.
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

t_sat = 373.16;
t_ratio = t_sat./T;
rt_ratio = 1.0./t_ratio;
sl_pressure = 1013.246;

c1 = 7.90298;
c2 = 5.02808;
c3 = 1.3816e-7;
c4 = 11.344;
c5 = 8.1328e-3;
c6 = 3.49149;

tmp = ( -1.0 * c1 * ( t_ratio - 1.0 ) ) + ...
      ( c2 * log10( t_ratio ) ) - ...
      ( c3 * ( 10.0.^( c4 * ( 1.0 - rt_ratio ) ) - 1.0 ) ) + ...
      ( c5 * ( 10.0.^( -1.0 * c6 * ( t_ratio - 1.0 ) ) - 1.0 ) ) + ...
      log10( sl_pressure );
svp = 10.0.^tmp;



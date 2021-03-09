%     Compute saturation vapor pressure (es,in mb) over water or ice at
%     temperature tk (kelvins), using the Goff-Gratch formulation (List,1963).
%     input arguments: 
%          tk   = temperature (K)
%          rh   = relative humidity (fraction)
%          ice  = switch to calculate saturation vapor pressure over
%                 water only (0) or water and ice, depending on tk (1)
%     output arguments: 
%          e    = vapor pressure (mb)
%          rho  = vapor density (g/m3)

function [e,rho] = Vapor_xxx(tk,rh,ice);

rvap = constants('Rwatvap'); % [J kg-1 K-1]
rvap = rvap * 1e-5;          % [J kg-1 K-1] -> [hPa * m2 g-1 K-1]
  
% if ( (tk > 263.16) | (ice==0) )
%    % for water...
%    y = 373.16 ./ tk;
%    es = -7.90298 * (y-1.) + 5.02808 * log10(y) -...
%          1.3816e-7 * (10 .^ (11.344 * (1. - (1./ y))) - 1.) +...
%          8.1328e-3 * (10 .^ (-3.49149 * (y - 1.)) - 1.) +...
%          log10(1013.246);
% else
%    % for ice...
%    y = 273.16 ./ tk;
%    es = -9.09718 * (y - 1.) - 3.56654 * log10(y) +...
%          0.876793 * (1.- (1. ./ y)) + log10(6.1071);
% end
  
% over water...
y = 373.16 ./ tk;
es = -7.90298 * (y-1.) + 5.02808 * log10(y) -...
     1.3816e-7 * (10 .^ (11.344 * (1. - (1./ y))) - 1.) +...
     8.1328e-3 * (10 .^ (-3.49149 * (y - 1.)) - 1.) +...
     log10(1013.246);

if ice == 1
   % over ice if tk < 263.16
   indx = find( tk < 263.16 );
   y = 273.16 ./ tk(indx);
   es(indx) = -9.09718 * (y - 1.) - 3.56654 * log10(y) +...
               0.876793 * (1.- (1. ./ y)) + log10(6.1071);
end

 es = 10. .^ es;
  
%     Compute vapor pressure and vapor density.
%     The vapor density conversion follows the ideal gas law: 
%     vapor pressure = vapor density * rvapor * tk
  
e = rh .* es;
rho = e ./ (rvap * tk);


return
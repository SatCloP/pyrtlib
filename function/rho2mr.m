%
% function mr = rho2mr(rho,t,p);
%
% determine water vapor mass mixing ratio (g/kg) given
% reference pressure (mbar), temperature (t,K), and
% water vapor density (g/m3).
%
% Equations were provided by Holger Linne' from Max Planck Institute.
%
% Nico 2002/05/09 (Looking at rho2mr.email.m from DCT)
% Nico 2018/09/26

function [mr] = rho2mr(rho,tk,p)

   mr = rho .* (tk) ./ (p * 0.3477);

   % Nico 2018/09/26
   % I think the above is an approximation valid within ~1%.
   % See mr2rho.m
   rvap = constants('Rwatvap') * 1e-5; % [J kg-1 K-1] -> [hPa * m2 g-1 K-1]
   eps = 0.621970585;                  % Rdry/Rvap
   mr = rho .* (1e3*eps*rvap) ./ (p./tk - rvap * rho);
   % where mr is the water vapor mass mixing ratio, i.e. Mv/Md (not the specific humidity, q=Mv/Mt=(mr/(1+mr)))
      
return


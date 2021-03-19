%
% function [rho] = mr2rho(mr,t,p);
%
% determine water vapor density (g/m3) given
% reference pressure (mbar), temperature (t,K), and
% water vapor mass mixing ratio (g/kg)
%
% Equations were provided by Holger Linne' from Max Planck Institute.
%
% Nico 2002/05/09 (Looking at rho2mr.email.m from DCT)
% Nico 2018/06/20

function [rho] = mr2rho(mr,tk,p);

     rho = mr .* p .* 0.3477 ./ (tk);

     % Nico 2018/06/20
     % I think the above is an approximation valid within ~1%.
     % To be consistent with Vapor_xxx.m and mr2rh.m (see
     % Compute_Transmittances_for_RTTOV_dsb.m), it should be:
     rvap = constants('Rwatvap') * 1e-5; % [J kg-1 K-1] -> [hPa * m2 g-1 K-1]
     eps = 0.621970585;                  % Rdry/Rvap
     rho = mr .* p .* 1./((1e3*eps+mr)*rvap) ./ (tk);
     % where mr is the water vapor mass mixing ratio, i.e. Mv/Md (not the specific humidity, q=Mv/Mt=(mr/(1+mr)))
     
return

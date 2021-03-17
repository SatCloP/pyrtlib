% This function computes H2O partial pressure (e,mbar) given pressure (p,mbar), temperature (t,K), and water vapor density (g/m3)
% It calls rho2mr and mr2e.
%
% NB: An approximation (not used here) is e = rho * t / 217 (this is used by Rosenkranz in his model)
%
%  USAGE:
%  e = rho2e(rho,t,p)
% 
%  History:
%      2015/12/14 - first version (Nico)

function e = rho2e(rho,t,p)

% computes H2O mass mixing ratio (g/kg) given pressure (p,mbar), temperature (t,K), and water vapor density (g/m3).
mr = rho2mr(rho,t,p);

% computes H2O partial pressure (e,mbar) given pressure (p,mbar) and H2O mass mixing ratio (mr,g/kg)
e = mr2e(p,mr);

return
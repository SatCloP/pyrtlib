function e = mr2e(p,mr);

%
% function e = mr2e(p,mr);
%
% compute H2O partial pressure (e,mbar) given
% pressure (p,mbar) and H2O mass mixing ratio (mr,g/kg)
%
% DCT 3/6/00
%

% ratio of water mass to dry air mass
eps = 0.621970585;

e = p.*mr./(1000.*eps + mr);

function mr = e2mr(p,e);

%
% function mr = e2mr(p,e);
%
% compute H2O mass mixing ratio (mr,g/kg) given
% pressure (p,mbar) and H2O partial pressure (e,mbar)
%
% DCT 3/6/00
%

% ratio of water mass to dry air mass
eps = 0.621970585;

mr = 1000.*eps.*e./(p-e);


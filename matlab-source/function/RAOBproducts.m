% This function computes derivate products
% from a standerd radiosonde output (z,p,t,rh)
% [e,tv,pt,N,vpt] = RAOBproducts(z,p,t,rh); km,mb,k,%/100
%
% History
% 2012/04/03 - Modified to include refractive index N

function [e,tv,pt,N,vpt] = RAOBproducts(z,p,t,rh)

addpath('C:/Matools/FROMDCT/');

% Constants
eps = 0.622;
p0 = 1000; % reference pressure
a = 0.286;

% Compute partial pressure
e = rh2e(p,t,rh);

% Compute virtual temperature
dum = (e./p)*(1-eps);
tv = t ./ (1-dum); 

% Compute potential temperature
pt = t .* (p0./p).^a;

% Compute virtual potential temperature
vpt = pt ./ (1-dum); 

% Compute refractive index
% Reference: Leroy et al, Testing climate models using GPS radio occultation: A sensitivity analysis, JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 111, D17105, doi:10.1029/2005JD006145, 2006
a = 77.6; % K/hPa
b = 373*1e3; % K^2/hPa
N = a*( p ./ t) + b * e ./ (t.^2);

return
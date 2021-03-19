% Quasi-gaussian distribution with up-to-third higher order terms.  
% Reference: "Scanning-laser glint measurements of sea-surface slope statisics" 
%             Shaw and Churnside, Applied Optics, V.36, n.18, 1997
%
% function [PD]=gaussdho(x,x0,sgm,skewness,kurtosis)

function [PD]=gaussdho(x,x0,sgm,c3,c4)

c0 = 1;

eta = (x-x0) / sgm ;

%normc = 1 / ( sqrt(2*pi) );
normc = 1 / ( sqrt(2*pi*sgm^2) ); % altrimenti non e' ortonormale
exarg = (eta).^2 / 2;

Gx = normc * exp(-exarg);

H0 = 1;
H3 = eta.^3 - 3*eta;
H4 = eta.^4 - 6*(eta.^2) + 3;

A0 = c0 * H0 .* Gx;
A3 = (1/6) * c3 * H3 .* Gx;
A4 = (1/24) * c4 * H4 .* Gx;

PD = A0 + A3 + A4;

return
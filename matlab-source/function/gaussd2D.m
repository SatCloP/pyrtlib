% FUNCTION GAUSSD2D
%
% Two dimensional Gaussian distribution.
% It takes in input two uncorrelated variables (the first can be a vector,
% while the second has to be a scalar) with their own center of distribution 
% and standard deviation.
% It gives in output the Multivariate Normalized Gaussian Distribution.
%
% Ref: "Multivariate Statistical Methods", D.F. Morrison, 
% McGraw-Hill series in Probability and Statistics, pagg. 80/81
%
% Es:
%    [gd]=gaussd2D(x,y,cx,cy,sx,sy)
%    dblquad('gaussd2D',-0.3,0.3,-0.3,0.3,0.001,'quad8')
%
% Inputs:
%         x: first dimension variable (can be a vector)
%         y: second dimension variable (must be a scalar)
%        cx: first dimension central value (defualt 0)
%        cy: second dimension central value (defualt 0)
%        sx: first dimension standard deviation (defualt 1)
%        sy: second dimension standard deviation (defualt 1)
% Outputs:
%        gd: Two dimensional Normalized Gaussian Distribution (1,length(x))
%
% Nico, Jun 2001

function [gd]=gaussd2D(x,y,cx,cy,sx,sy)

% Questo serve solo per la prova della normalizzazione
if nargin == 2   
   cx=0; cy=0;
   sx=0.05; sy=0.03;
end

norm = (2*pi) * sx * sy;
norm = 1 / norm;  

arg  = 0.5 * ( ((x-cx)/sx).^2 + ((y-cy)/sy)^2 ) ; 
gd = norm * exp(-arg);

return


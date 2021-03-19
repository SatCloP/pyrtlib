% FUNCTION GAUSSD
%
% It takes in input the vector of x, the center of distribution 
% and the standard deviation with respect to the center. 
% It gives in output the Normalized Gaussian Distribution over x.
%
% N.B: -3dB full-width gives a good approximation of 2*sigma, but IT'S NOT 2*sigma.
%      In a gaussian curve, sigma is the distance between simmetry axes and 
%      the flex-point. According to this, -3db corrispondes to 2.3548*sigma, 
%      as comes out solving for x the equation exp(-(x/2sigma)^2)=1/2
%
%      Matlab statistics toolbox has normpdf(x,x0,sigma)
%
% function [gd]=gaussd(x,c,s)
%
% Nico,1999

function [gd]=gaussd(v,c,s)

  norm = ( (2*pi)^(1/2) ) * s;
  norm = 1 / norm;
  
  arg  = ( ( (v-c)/s ).^2 ) / 2;
  
  gd   = norm * exp(-arg);

return
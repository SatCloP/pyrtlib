% THIS FUNCTION GIVES BACK BRIGHTNESS TEMPERATURE ONCE
% GIVEN IN INPUT IR RADIANCE AND WAVE NUMBER.
% THE INPUT UNITS SHOULD BE [W/(m2 sr cm-1)] AND [cm-1]
% THE OUTPUT UNITS WILL BE [K]
%
% Inputs:
%         rad = Radiance [W/(m2 sr cm-1)]
%         wn  = Wavenumber [cm-1]
% Outputs:
%         Tb  = Brightness Temperature [K] 
% Es:
%    [Tb]=radian2tb(rad,wn);

function [Tb]=radian2tb(rad,wn)

  h=constants('planck');    % [J Hz-1] Planck constant 
  K=constants('boltzmann'); % [J K-1] Boltzmann constant
  c=constants('light');     % [m s-1] Light speed
  a=1/100;                  % Conversion factor

  ONE=(h*c*wn)/(a*K);
  ADD=(2*h*(c^2)*wn.^3)./(a^4.*rad);
  TWO=log(1+ADD);
  
  Tb=ONE./TWO;

return
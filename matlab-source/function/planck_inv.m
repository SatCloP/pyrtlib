% THIS FUNCTION COMPUTES BRIGHTNESS TEMPERATURE PROVIDING IN INPUT 
% RADIANCE AND THE FREQUENCY (OR WAVENUMBER OR WAVELENGTH) VECTOR.
%
% Usage:
%          Tb = planck_inv(rad,x,units);
% Inputs:
%         rad: radiance [W/(m2 sr cm-1) or W/(m2 sr m) or W/(m2 sr Hz) or K]
%           x: wavenumber [cm-1] or wavelength [m] or frequency [Hz]
%       units: string  
%              'k' wave number (cm-1)  input radiance in [ W/(m2 cm-1 sr) ]
%              'l' wavelength (m)      input radiance in [ W/(m2 m sr) ]
%              'f' frequency (Hz)      input radiance in [ W/(m2 Hz sr) ]
%              'RJE' frequency (Hz)    input radiance in [ K ] using Rayleigh-Jeans Equivalent Tb scaling (Janssen, 1993, pag. 6)
% Outputs:
%          Tb: Thermodinamic Brightness Temperature [K] 
% Es:
%          Tb = planck_inv(rad,freq,'f');
%
% see also planck_f, planck_approx
%
% Nico Nov, 2005
%
% History:
% 2005/11/02 - first version

function Tb = planck_inv(rad,x,units)

h = constants('planck');    % [J Hz-1] Planck constant 
K = constants('boltzmann'); % [J K-1] Boltzmann constant
c = constants('light');     % [m s-1] Light speed
a = 1/100;                  % Conversion factor

switch units
   case 'k'   
     % INPUT UNITS: W/(m2 cm-1 sr)
     first = (h*c*x)/(a*K);
     addnd = (2*h*(c^2)*x.^3)./(a^4.*rad);
     second = log(1+addnd);  
   case 'l'   
     % INPUT UNITS: W/(m2 m sr)
     fprintf(1,'\n ATTENZIONE!: QUESTA FUNZIONE POTREBBE DARE RISULTATI \n NON ACCURATI A CAUSA DELLA PRECISIONE DI CALCOLO!\n')
     first = h*c./(K*x);
     addnd = (2*h*(c^2))./(rad.*(x.^5));
     second = log(1+addnd);  
   case 'f'
     % INPUT UNITS: W/(m2 Hz sr)
     first = h*x/K;
     addnd = (2*h*x.^3/(c^2))./(rad);
     second = log(1+addnd); 
   case 'RJE'   
     % INPUT UNITS: K
     first = h*x/K;
     addnd = (h*x)./(K*rad);
     second = log(1+addnd);  
end

Tb = first./second;
  
return
% FUNCTION THAT CALCULATES THE PLANCK FUNCTION, 
% PROVIDING IN INPUT THE PHYSICAL TEMPERATURE (K) AND
% THE FREQUENCY (OR WAVENUMBER OR WAVELENGTH) VECTOR.
% YOU NEED TO SPECIFY INPUT UNITS, CHOOSING BETWEEN:
%
% 'k' wave number (cm-1)  output in [ W/(m2 cm-1 sr) ]
% 'l' wavelength (m)      output in [ W/(m2 m sr) ]
% 'f' frequency (Hz)      output in [ W/(m2 Hz sr) ]
%
% ES: [B_f]=planck_f(T,freq,'f')

function [B]=planck_f(T,x,units)

h=constants('planck');    % [J Hz-1] Planck constant 
K=constants('boltzmann'); % [J K-1] Boltzmann constant
c=constants('light');     % [m s-1] Light speed
a=1/100; % Conversion factor

switch units
   
case 'k'
   
   % OUTPUT UNITS: W/(m2 cm-1 sr)
   esp=(h*c*x)./(a*K*T);
   first=2*h*(c^2)*(x.^3)/(a^4);
   second=(1./(exp(esp)-1));
   
case 'l'
   
   % OUTPUT UNITS: W/(m2 m sr)
   fprintf(1,'\n ATTENZIONE!: QUESTA FUNZIONE POTREBBE DARE RISULTATI \n NON ACCURATI A CAUSA DELLA PRECISIONE DI CALCOLO!\n')
   esp=(h*c)./(K*T.*x);
   first=(2*h*c^2)./(x.^5);
   second=(1./(exp(esp)-1));

case 'f'
   
   % OUTPUT UNITS: W/(m2 Hz sr)
   esp=(h*x)./(K*T);
   first=2*h*(x.^3)/(c^2);
   second=(1./(exp(esp)-1));
   
end

B=first.*second;


return
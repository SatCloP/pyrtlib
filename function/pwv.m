%
% function pw = pwv(pres,w)
%
% Computes the integrated water vapor content (or precipitable water vapor PWV)
% (cm) giving in input the atmospheric profiles of pressure levels (mb) 
% and water vapor mass mixing ratio (g/Kg).
%
% DCT

function pw = pwv(pres,w)

   g = 9.80665;     % gravity
   pres = flipud(pres)*100; % Pa
   w = flipud(w)/1000;      % g/g or kg/kg
   sp_hum = w./(1+w);
   %pw = (0.1/g)*integral(pres,sp_hum);
   pw = (0.1/g)*trapz(pres,sp_hum);

return
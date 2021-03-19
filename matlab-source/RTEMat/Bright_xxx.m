% ...................................................................
%     Function to compute temperature from the modified Planck
%     radiance (Planck function without the constants 2h(v^3)/(c^2).
%
%     inputs passed as arguments: 
%
%        hvk  =   [Planck constant (J*S)] * [frequency (Hz)]
%		  --------------------------------------
%                     [Boltzmann constant (J/K)]
%        boft =   modified Planck radiance -
%		  (equation (4) from Schroeder & Westwater, 1991)
% ...................................................................
function Tb = Bright_xxx(hvk,boft);

Tb = hvk ./ log(1.0 + (1.0 ./ boft));

return

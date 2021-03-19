% This script computes the vertical gradient of the potential temperature,
% also called dtheta/dz, providing in input vertical profiles of P [mb], 
% T [K], Z [m]
%
% USAGE:
% [dthdz,zhp,theta] = dthetadz(Pmb,TK,Zm);
% 
% INPUT:
% vertical profiles of 
% P [mb], T [K], Z [m]
%
% OUTPUT:
% vertical profiles of 
% dthdz [K/m], defined on zhp heights
% zhp [m] heights of half-point (to compute gradients)
% theta [K], defined on Zm

function [dthdz,zhp,theta] = dthetadz(Pmb,TK,Zm)

    es = 0.286; % Exponent
    %fact1 = (Pmb(1)*ones(size(Pmb))./Pmb).^es; 
    fact1 = (1000*ones(size(Pmb))./Pmb).^es;  % See paper from B. B. Stankov et al., 2002, JTECH (ma non cambia molto...)
    theta = TK .* fact1; % Potential Temperature

    [dthdz,zhp] = vgrad(Zm,theta);

return


% vgrad compute vertical gradients
function [dxdz,zm] = vgrad(z,x)

N = length(z);

for i = 1:N-1
    
    dxdz(i) = ( x(i+1) - x(i) ) / ( z(i+1) - z(i) );
    zm(i) = ( z(i+1) + z(i) ) / 2;

end

return
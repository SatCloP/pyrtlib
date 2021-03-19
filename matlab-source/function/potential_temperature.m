% This function computes the potential temperature theta 
% from vertical profiles of P [mb], T [K], Z [m]
%
% USAGE:
%         theta = potential_temperature(Pmb,TK,Zm);
% 
% INPUT: 
%         P [mb] - vertical profile
%         T [K] - vertical profile
%         Z [m] - vertical profile
%
% OUTPUT:
%         theta [K] - vertical profile of theta defined on Zm
%
% 2021/02/26 - Nico: first version, looking at dthetadz.m

function theta = potential_temperature(Pmb,TK,Zm)

    es = 0.286; % Exponent
    %fact1 = (Pmb(1)*ones(size(Pmb))./Pmb).^es; 
    fact1 = (1000*ones(size(Pmb))./Pmb).^es;  % See paper from B. B. Stankov et al., 2002, JTECH (ma non cambia molto...)
    theta = TK .* fact1; % Potential Temperature

return

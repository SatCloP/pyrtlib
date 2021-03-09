%   Computes the modified planck function (equation (4) in schroeder and
%   westwater, 1992: guide to passive microwave weighting function
%   calculations) for the cosmic background temperature, the mean radiating
%   temperature, and a profile of the atmospheric integral with and without
%   the cosmic background. Also computes an integral profile of atmospheric
%   absorption. For the integral profiles, the value at profile level i
%   represents the integral from the antenna to level i.
%   Also returns the cosmic background term for the rte.
%
%   inputs passed as arguments: 
%        frq     = channel frequency (GHz)
%        nl      = number of profile levels
%        tk      = temperature profile (K)
%        taulay  = profile of absorption integrated over each layer (np)
%   outputs: 
%        hvk     = [planck constant * frequency] / boltzmann constant
%        boft    = modified planck function for raob temperature profile
%        bakgrnd = background term of radiative transfer equation
%        boftatm = array of atmospheric planck radiance integrated (0,i)
%        boftotl = total planck radiance from the atmosphere plus bakgrnd
%        boftmr  = modified planck function for mean radiating temperature
%        tauprof = array of integrated absorption (np; 0,i)
function  [boftotl_sat,boftatm_sat,boftmr_sat,tauprof_sat,hvk,boft_sat,bakgrnd_sat] = Planck_xxx(frq,tk,taulay,Es);

Tc = constants('Tcosmicbkg');
Ts = tk(1); % Tsurface is ~ temperature of lowermost layer
h = constants('planck');
k = constants('boltzmann');
fHz = frq*1e9; % GHz -> Hz
hvk = fHz * h / k;
% maximum absolute value for exponential function argument
expmax = 125.0;
nl = length(tk);
tauprof = zeros(size(taulay));
boftatm = zeros(size(taulay));

tauprof_sat = zeros(size(taulay));
boftatm_sat = zeros(size(taulay));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First compute downwelling radiance (to be reflected by the surface)
% as in original Planck_xxx.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Es == 1;
   % no need to compute downwelling radiance
   boftotl = 0.0;
else
    boft(1) = TK2B_mod(hvk,tk(1));
    for i = 2:nl
        boft(i) = TK2B_mod(hvk,tk(i));
        boftlay = (boft(i-1) + boft(i) * exp(-taulay(i))) / (1.0 + exp(-taulay(i)));
        batmlay = boftlay * exp(-tauprof(i-1)) * (1.0 - exp(-taulay(i)));
        boftatm(i) =  boftatm(i-1) + batmlay;
        tauprof(i) = tauprof(i-1) + taulay(i);
    end
    %     compute the cosmic background term of the rte; compute total planck
    %     radiance for atmosphere and cosmic background; if absorption too large
    %     to exponentiate, assume cosmic background was completely attenuated.
    if tauprof(nl) < expmax
        boftbg = TK2B_mod(hvk,Tc);
        bakgrnd = boftbg * exp(-tauprof(nl));
        boftotl = bakgrnd + boftatm(nl);
        boftmr = boftatm(nl) / (1. - exp (-tauprof(nl)));
    else
        bakgrnd = 0.;
        boftotl = boftatm(nl);
        boftmr  = boftatm(nl);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Then compute upwelling radiance 
% Adapted from Planck_xxx.m, but from Satellite i-1 becomes i+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boft_sat(nl) = TK2B_mod(hvk,tk(nl));
for i = nl-1 : -1 : 1
    boft_sat(i) = TK2B_mod(hvk,tk(i));
    boftlay_sat = (boft_sat(i+1) + boft_sat(i) * exp(-taulay(i))) / (1.0 + exp(-taulay(i)));
    batmlay_sat = boftlay_sat * exp(-tauprof_sat(i+1)) * (1.0 - exp(-taulay(i)));
    boftatm_sat(i) = boftatm_sat(i+1) + batmlay_sat;
    tauprof_sat(i) = tauprof_sat(i+1) + taulay(i);
end

% The background is a combination of surface emission and downwelling 
% radiance (boftotl) reflected by the surface
if tauprof_sat(1) < expmax
    boftbg_sat  = Es * TK2B_mod(hvk,Ts) + (1-Es) * boftotl;      % SAT: eps * B(Tsrf) + (1-eps) B_dw
    %boftbg_sat  = Es * TK2B_mod(hvk,Ts);      % SAT: eps * B(Tsrf) + (1-eps) B_dw
    bakgrnd_sat = boftbg_sat * exp(-tauprof_sat(1));                 % SAT: nl -> 1
    boftotl_sat = bakgrnd_sat + boftatm_sat(1);                  % SAT: nl -> 1
    boftmr_sat  = boftatm_sat(1) / (1. - exp (-tauprof_sat(1))); % SAT: nl -> 1
else
    bakgrnd_sat = 0.;
    boftotl_sat = boftatm_sat(1);
    boftmr_sat  = boftatm_sat(1);
end

return



function Btilde = TK2B_mod(hvk,T) 

Btilde = 1.0 / (exp(hvk/T) - 1.0);

return

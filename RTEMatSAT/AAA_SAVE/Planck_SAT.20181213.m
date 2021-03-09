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
function  [boftotl,boftatm,boftmr,tauprof,hvk,boft,bakgrnd] = Planck_xxx(frq,tk,taulay,Es);

%Tbck = constants('Tcosmicbkg');
Tbck = tk(1); % background now is ~ surface temperature
h = constants('planck');
k = constants('boltzmann');
fHz = frq*1e9; % GHz -> Hz
hvk = fHz * h / k;
% maximum absolute value for exponential function argument
expmax = 125.0;
nl = length(tk);
tauprof = zeros(size(taulay));
boftatm = zeros(size(taulay));

% boft(1) = TK2B_mod(hvk,tk(1));    
% for i = 2:nl
%     boft(i) = TK2B_mod(hvk,tk(i));
%     boftlay = (boft(i-1) + boft(i) * exp(-taulay(i))) / (1.0 + exp(-taulay(i)));
%     batmlay = boftlay * exp(-tauprof(i-1)) * (1.0 - exp(-taulay(i)));
%     boftatm(i) =  boftatm(i-1) + batmlay;
%     tauprof(i) = tauprof(i-1) + taulay(i);
% end
%
% From Satellite i-1 becomes i+1
boft(nl) = TK2B_mod(hvk,tk(nl));
for i = nl-1 : -1 : 1
    boft(i) = TK2B_mod(hvk,tk(i));
    boftlay = (boft(i+1) + boft(i) * exp(-taulay(i))) / (1.0 + exp(-taulay(i)));
    batmlay = boftlay * exp(-tauprof(i+1)) * (1.0 - exp(-taulay(i)));
    boftatm(i) =  boftatm(i+1) + batmlay;
    tauprof(i) = tauprof(i+1) + taulay(i);
end

%     compute the cosmic background term of the rte; compute total planck
%     radiance for atmosphere and cosmic background; if absorption too large
%     to exponentiate, assume cosmic background was completely attenuated.
% if tauprof(nl) < expmax
%     boftbg = TK2B_mod(hvk,Tbck);
%     bakgrnd = boftbg * exp(-tauprof(nl));
%     boftotl = bakgrnd + boftatm(nl);
%     boftmr = boftatm(nl) / (1. - exp (-tauprof(nl)));
% else
%     bakgrnd = 0.;
%     boftotl = boftatm(nl);
%     boftmr  = boftatm(nl);
% end
%
% From Satellite i-1 becomes i+1
%     compute the surface background term of the rte; compute total planck
%     radiance for atmosphere and surface background; if absorption too large
%     to exponentiate, assume surface background was completely attenuated.
if tauprof(nl) < expmax
    boftbg = TK2B_mod(hvk,Tbck);
    bakgrnd = Es * boftbg * exp(-tauprof(1));       % SAT: nl -> 1 * Surface emissivity!
    boftotl = bakgrnd + boftatm(1);                 % SAT: nl -> 1
    boftmr = boftatm(1) / (1. - exp (-tauprof(1))); % SAT: nl -> 1
else
    bakgrnd = 0.;
    boftotl = boftatm(1);
    boftmr  = boftatm(1);
end

return

function Btilde = TK2B_mod(hvk,T) 

Btilde = 1.0 / (exp(hvk/T) - 1.0);

return

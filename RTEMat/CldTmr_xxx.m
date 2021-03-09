% ...................................................................
%     Computes the mean radiating temperature of a cloud with base and top at
%     profile levels ibase and itop, respectively.  The algorithm assumes that
%     the input cloud is the lowest (or only) cloud layer observed.
%     If absorption is not too big, compute tmr of lowest cloud layer (base
%     at level ibase, top at level itop). Otherwise, set error flag and return. 
% *** NOTE: This algorithm is not designed for multiple cloud layers.***
%
%     inputs passed as arguments: 
%          ibase   = profile level at base of lowest cloud
%          itop    = profile level at top of lowest cloud
%          hvk     = (planck constant * frequency) / boltzmann constant
%          tauprof = integral profile of absorption (np; i = integral (0,i))
%          boftatm = integral profile of atmospheric planck radiance
% *** note: hvk, tauprof, and boftatm can be obtained from subroutine planck.
%
%     output: 
%          tmrcld = tmr of lowest cloud layer (k) 
% ...................................................................
function tmrcld = CldTmr_xxx(ibase,itop,hvk,tauprof,boftatm);

% maximum absolute value for exponential function argument
expmax = 125.0;

% ... check if absorption too large to exponentiate...
if tauprof(ibase) > expmax
   display('from CldTmr_xxx: absorption too large to exponentiate for tmr of lowest cloud layer');
   return
end

% compute radiance (batmcld) and absorption (taucld) for cloud layer.
% (if taucld is too large to exponentiate, treat it as infinity.) 
batmcld = boftatm(itop) - boftatm(ibase);
taucld = tauprof(itop) - tauprof(ibase);
if taucld > expmax
    boftcld = batmcld * exp(tauprof(ibase));
else
    boftcld = batmcld * exp(tauprof(ibase)) / (1. - exp (-taucld));
end
	
% compute cloud mean radiating temperature (tmrcld)
tmrcld = Bright_xxx(hvk,boftcld);

return
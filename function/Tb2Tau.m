% Function Tb2Tau converts Brightness Temperature to Tau or inverse, 
% depending on inputs.
% If the frequency is provided in input, the exact formulation or the
% approximation suggested in Janssen (eq 1.21), is computed.
%
% Usage:
%        Out = Tb2Tau(In,Tmr,flag,frq,jan);
% Inputs:
%        In: Tb (K) or Tau (Np)
%       Tmr: Mean Radiative Temperature (K)
%      flag: 1 for Tb->Tau, 0 for Tau->Tb 
% Optional input:
%       frq: frequency (GHz). If this is given, the exact formulation
%            using the Planck function, is computed
%       jan: if 1, computes tau using Janssen approximation 
%            tau = log((Tmr-Tbc)/(Tmr-Tb)) where
%            Tbc = (hf/2k)(exp(hf/kTbg)+1)/(exp(hf/kTbg)-1) and Tbg = 2.736
% Output:
%       Out: Tau (Np) or Tb (K)
% Es:
%     Tau = Tb2Tau(Tb,Tmr,1);
%     Tau = Tb2Tau(Tb,Tmr,1,frq);
%     Tau = Tb2Tau(Tb,Tmr,1,frq,1);
%     Tb = Tb2Tau(Tau,Tmr,0);
%     Tb = Tb2Tau(Tau,Tmr,0,frq);
%     Tb = Tb2Tau(Tau,Tmr,0,frq,1);
%
% Nico, ??
%
% History:
% 2005/11/02 - Included the choice for exact formulation
% 2005/11/04 - Included the choice for Janssen approximation

function [Out,compute,Bv_tbt,Bv_tmr,Bv_tbg,Tcos] = Tb2Tau(In,Tmr,flag,frq,jan)

Tcos = constants('Tcosmicbkg'); % 2.736 +/- 0.017 K (from Janssen, Atmospheric Remote Sensing by Microwave Radiometry, pag.12)
%Tcos = 2.75; % I used this for L24 data

if nargin < 4
   compute = 'apprx';
   Bv_tbt=[]; Bv_tmr=[]; Bv_tbg=[];
elseif nargin == 4
   compute = 'exact';
   frqHz = 1e9*frq';   
   Bv_tmr = planck_f(Tmr,frqHz,'f');
   Bv_tbg = planck_f(Tcos,frqHz,'f');
else
   compute = 'janss';
   h = constants('planck');    % [J Hz-1] Planck constant 
   k = constants('boltzmann'); % [J K-1] Boltzmann constant
   c = constants('light');     % [m s-1] Light speed
   frqHz = 1e9*frq';      
   expn = exp( h*frqHz./(k*Tcos) );
   Tcos = h*frqHz./(2*k) .* ( (expn+1)./(expn-1) );
   Bv_tbt=[]; Bv_tmr=[]; Bv_tbg=[];
end


if flag   
     
   Tb = In;
   switch compute
       case 'apprx'
            Tau = log ((Tmr-Tcos)./(Tmr-Tb));
       case 'exact'
            Bv_tbt = planck_f(Tb,frqHz,'f');
            Tau = log( (Bv_tmr - Bv_tbg) ./ (Bv_tmr - Bv_tbt) );
       case 'janss'
            Tau = log ((Tmr-Tcos)./(Tmr-Tb)); % Tcos depends on frequency now
   end
   if 1 % If 1, the following avoids complex numbers, but sometimes it's better to know when it's going crazy
      jx = find( Tmr-Tb <= 0 );
      Tau(jx) = NaN;
   end
   Out = Tau;
    
 else
     
   Tau = In;
   switch compute
       case 'apprx'
            Tb = Tmr - ( Tmr-Tcos ) .* exp(-Tau);
       case 'exact'
            BTb = Bv_tmr - (Bv_tmr - Bv_tbg) .* exp(-Tau);
            Tb = planck_inv(BTb,frqHz,'f');
       case 'janss'
            Tb = Tmr - ( Tmr-Tcos ) .* exp(-Tau); % Tcos depends on frequency now
   end
   Out = Tb;
   
end 
   
return
% ...................................................................
%     Multiplies cloud density profiles by a given fraction and computes the
%     corresponding cloud liquid and ice absorption profiles, using Rosenkranz's
%     cloud liquid absorption routine ABLIQ and ice absorption of Westwater
%     [1972: Microwave Emission from Clouds,13-14].   ---- Yong Han, 4/20/2000
%
%     inputs passed as arguments: 
%          tk      = temperature profile (k)
%          denl    = liquid density profile (g/m3) 
%          deni    = ice density profile (g/m3)
%          frq     = frequency array (GHz)
%     outputs: 
%          aliq  = liquid absorption profile (np/km)
%          aice  = ice absorption profile (np/km) 
% ...................................................................

function [aliq,aice] = CldAbs_xxx(tk,denl,deni,frq);


nl = length(tk);
c = 2.99792458e10; % cm/s
ghz2hz = 1e9;
db2np = log (10.) * 0.1; % conversion factor: decibels to nepers
wave = c / (frq * ghz2hz); % cm/2 * Hz = cm

aliq = zeros(size(denl));
aice = zeros(size(deni));

for i = 1:nl
  
    % Compute liquid absorption np/km.
    if denl(i) > 0
       aliq(i) = ABLIQ(denl(i),frq,tk(i));
    end
  
    % compute ice absorption (db/km); convert non-zero value to np/km.
    if deni(i) > 0
       aice(i) = (8.18645 / wave) * deni(i) * 9.59553e-4;
       aice(i) = aice(i) * db2np;
    end
    

end

return


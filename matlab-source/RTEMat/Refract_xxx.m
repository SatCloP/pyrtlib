% ...................................................................
%     Computes profiles of wet refractivity, dry refractivity, 
%     refractive index.  Refractivity equations were taken from G.D. 
%     Thayer, 1974:  An improved equation for the radio refractive 
%     index of air. Radio Science, vol.9,no.10, 803-807.
% *** These equations were intended for frequencies under 20 GHz ***
%
%     inputs: 
%          p   = pressure profile (mb)
%          tk  = temperature profile (K)
%          e   = vapor pressure profile (mb)
%     outputs: 
%          dryn    = dry refractivity profile
%          wetn    = wet refractivity profile
%          refindx = refractive index profile
%    Example:
%          [dryn,wetn,refindx] = Refract_xxx(p,tk,e);
% ...................................................................

function [dryn,wetn,refindx] = Refract_xxx(p,tk,e)

nl = length(p);

for i = 1:nl
    
    % Calculate dry air pressure (pa) and celsius temperature (tc).    
    pa = p(i) - e(i);
    tc = tk(i) - 273.16;
    
    % Calculate inverse wet (rzw) and dry (rza) compressibility factors.
    
    tk2 = tk(i) * tk(i);
    tc2 = tc * tc;
    rza = 1. + pa * (57.90e-8 * (1.+0.52/tk(i)) - 9.4611e-4 * tc/tk2);
    rzw = 1. + 1650. * (e(i) / (tk(i)*tk2)) * (1.-0.01317 * tc + 1.75e-4 * tc2 + 1.44e-6 * (tc2 * tc));
    
    % Calculate wet refractivity, dry refractivity, and refractive index.
    
    wetn(i) = (64.79 * (e(i)/tk(i)) + (3.776e+5) * (e(i)/tk2)) * rzw;
    dryn(i) = 77.6036 * (pa / tk(i)) * rza;
    refindx(i) = 1. + (dryn(i) + wetn(i)) * 1.e-6;
    
end


return

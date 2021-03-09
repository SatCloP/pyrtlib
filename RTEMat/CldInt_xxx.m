% ...................................................................
%     Integrates cloud water density over path ds (linear algorithm). 
%
%     inputs: 
%          dencld = cloud cloud water density profile (g/m3)
%          ds     = vector containing layer depth profiles (km)
%          nlay   = number of cloud layers in the profile
%          lbase  = array containing profile levels corresponding to cloud bases
%          ltop   = array containing profile levels corresponding to cloud tops 
%     outputs: 
%          scld   = integrated cloud water density (cm)
% ...................................................................

function scld = CldInt_xxx(dencld,ds,lbase,ltop);

ncld = length(lbase);

scld = 0.0;
for l = 1:ncld
    for i = lbase(l)+1:ltop(l)
        scld = scld + ds(i) * (0.5 * (dencld(i) + dencld(i-1)));
    end
end

% convert the integrated value to cm.
scld = scld * 0.1;

return
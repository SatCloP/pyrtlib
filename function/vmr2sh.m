% vmr2sh convert volume mixing ratio (ppm) in specific humidity (g/g)
%
% USAGE:
%      sh = vmr2sh(vmr);
%
% History:
%     2012/05/18 - first version (Nico)

function sh = vmr2sh(vmr)

  Md = 28.966; % molecular mass of dry air
  Mw = 18.016; % molecular mass of water

  sh = vmr*Mw ./ (vmr*Mw + 1e6*Md);
  
return

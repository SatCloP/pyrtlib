% sh2vmr converts specific humidity (g/g) in volume mixing ratio (ppm)
%
% USAGE:
%      vmr = sh2vmr(sh);
%
% History:
%     2012/05/18 - first version (Nico)

function vmr = sh2vmr(sh)

  Md = 28.966; % molecular mass of dry air
  Mw = 18.016; % molecular mass of water

  % complete formula (reference: Stefan Buhler presentation)
  vmr = 1e6 * sh ./ ( (1-sh)*Mw/Md + sh );

  % approximated formula 
  vmr0 = 1e6 * Md/Mw * (sh ./ (1-sh));
  % percentage of approximation
  approx = (vmr0-vmr)/vmr * 100;
  
return

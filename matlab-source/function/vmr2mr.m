% This function convert volume mixing ratio (ppmv) into mass mixing ratio (g/kg).
%
% Example:
%     mr = vmr2mr(vmr);
% 
% INPUT:
%     vmr: volume mixing ratio [ppm] e.g. V(H2O) / V(dry air)
% 
% OUTPUT:
%      mr: mass mixing ratio [g/kg] e.g. g(H2O) / kg(dry air)
% 
% REFERENCE:
%     mr = m_H2O/m_d = n_w*Mw / n_d*Md (n_w and n_d are the number of moles of water vapor and dry air, respectively)
%    vmr = n_w / n_d
%     mr = Mw/Md * vmr
%
% HISTORY:
%  2014/09/04 - First version (Nico)

function  mr = vmr2mr(vmr)
   
   Md = 28.966; % molecular mass of dry air
   Mw = 18.016; % molecular mass of water
   
   mr = Mw/Md * vmr * 1e-3; % mg/kg -> g/kg
   
return
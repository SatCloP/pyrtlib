% This function convert mass mixing ratio (g H2O per kg dry air) into 
% volume mixing ratio in terms of ppmv.
%
% Example:
%     vmr = mr2vmr(mr);
% 
% INPUT:
%      mr: mass mixing ratio [g/kg] e.g. g(H2O) / kg(dry air)
% 
% OUTPUT:
%     vmr: volume mixing ratio [ppm] e.g. V(H2O) / V(dry air)
% 
% REFERENCE:
%     mr = m_H2O/m_d = n_w*Mw / n_d*Md (n_w and n_d are the number of moles of water vapor and dry air, respectively)
%    vmr = n_w / n_d
%    vmr = Md/Mw * mr
%
% HISTORY:
%  2014/09/04 - First version (Nico)

function  vmr = mr2vmr(mr)
   
   Md = 28.966; % molecular mass of dry air
   Mw = 18.016; % molecular mass of water
   
   vmr = Md/Mw * mr * 1e3; % ppm
   
return

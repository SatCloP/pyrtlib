% This function gives back absorption coefficient for liquid and ice water,
% as reported in the open litterature.
%
% Nico, Sep 2006
%
% Usage:
%   [wl,ac] = water_abscoef(flag);
% Input:
%   flag: 'liq' or 'ice' water
% Outputs:
%     wl: wavelength (micron)   
%     ac: absorption coefficient (cm-1) 
% Example:
% [wl_liq,ac_liq] = water_abscoef('liq');
% [wl_ice,ac_ice] = water_abscoef('ice');
% norm = 1e4; % (micron to cm-1) and (cm-1 to micron-1)
% plot(norm./wl_liq,ac_liq/norm,'b',norm./wl_ice,ac_ice/norm,'c'); xlim([400 1400]); ylim([0 0.5]); 
% xlabel('Wavenumber [cm^{-1}]'); ylabel('Abs Coef [\mum^{-1}]'); legend('liquid','ice');
%
% References:
% Hale G. M. and M. R. Querry, "Optical constants of water in the 200nm to
%      200µm wavelength region," Appl. Opt., 12, 555--563, (1973).
% Warren S. G. , "Optical constants of ice from the ultraviolet to the
%      microwave," Appl. Opt., 23, 1026--1225, (1984).
%
% History:
% 2006/09/12 - First version

function [wlength,abscoef] = water_abscoef(flag);

switch flag
    case 'liq'
       datfile = 'hale73.dat';
    case 'ice'
       datfile = 'warren95.dat';
end

MTX = load(datfile);
wlength = MTX(:,1)/1000; % nm -> micron
abscoef = MTX(:,2); % cm-1

return
% This function compute units conversions.
%
% Usage:
%
%           out = units_conversion(inp,FromTo)
% 
% Available:
%           FromTo                inp                 out
%           ------------------------------------------------
%           'GHzTOm'              f(GHz)              wl(m)
%           'mTOGHz'              wl(m)               f(GHz)
%           'GHzTOcm-1'           f(GHz)              wn(cm-1)
%           'cm-1TOGHz'           wn(cm-1)            f(GHz)
%           'GHz/kPaTOMHz/Torr'   gamma(GHz/kPa)      gamma(MHz/Torr)
%           'MHz/TorrTOGHz/kPa'   gamma(MHz/Torr)     gamma(GHz/kPa)
%           'cm-1/atmTOMHz/Torr'  I(cm-1/atm)         I(MHz/Torr)
%           'MHz/TorrTOcm-1/atm'  I(MHz/Torr)         I(cm-1/atm)
%
% Examples:
%           wn_inverse_cm = units_conversion(f_GHz,'GHzTOcm-1')
%           f_GHz = units_conversion(wn_inverse_cm,'cm-1TOGHz')
%
% NB: See also Unit_conversion.docx !!!

function [out,units] = units_conversion(inp,FromTo)

switch FromTo
    
    case 'GHz/kPaTOMHz/Torr'

         % inp gamma in GHz/kPa 
         % out gamma in MHz/Torr
         mb2torr = 0.750062; % mb2torr conversion factor; MHz/Torr -> MHz/mb (the conversion factor is 1/torr2mb, i.e. mb2torr)
         out = (inp * 1e2);   % GHz/kPa -> MHz/mb
         out = out / mb2torr; % MHz/mb -> MHz/Torr
         units = 'MHz/Torr';
        
    case 'MHz/TorrTOGHz/kPa'

         % inp gamma in MHz/Torr  
         % out gamma in GHz/kPa
         mb2torr = 0.750062; % mb2torr conversion factor; MHz/Torr -> MHz/mb (the conversion factor is 1/torr2mb, i.e. mb2torr)
         out = inp * mb2torr; % MHz/Torr -> MHz/mb
         out = out * 1e-2;   % MHz/mb -> GHz/kPa
         units = 'GHz/kPa';
        
    case 'GHzTOm'
         % inp frequency in GHz 
         % out wavelength in m
         [c,units]=constants('light');
         out = c / (inp * 1e9); 
         units = 'm';
        
    case 'mTOGHz'
         % inp wavelength in m 
         % out frequency in GHz
         [c,units]=constants('light');
         out = (c / inp) / 1e9; 
         units = 'GHz';
        
    case 'GHzTOcm-1'
         % inp frequency in GHz 
         % out wavenumber in cm-1
         [c,units]=constants('light');
         c = c * 1e2; %m/s -> cm/s
         out = (inp * 1e9) / c; 
         units = 'cm-1';
        
    case 'cm-1TOGHz'
         % inp wavenumber in cm-1
         % out frequency in GHz 
         [c,units]=constants('light');
         c = c * 1e2; %m/s -> cm/s
         out = (inp * c) * 1e-9; 
         units = 'GHz';

    case 'cm-1/atmTOMHz/Torr'

         % inp line intensity in cm-1/atm  
         % out line intensity in MHz/Torr
         mb2torr = 0.750062;
         atm2torr = 1013 * mb2torr; % 1013 mb * mb2torr 
         [c,units]=constants('light');
         c = c * 1e2; %m/s -> cm/s
         out = (inp * c) * 1e-6; % cm-1*atm-1 -> MHz/atm
         out = out /atm2torr;    % MHz/atm -> MHz/Torr
         units = 'MHz/Torr';

    case 'MHz/TorrTOcm-1/atm'

         % inp line intensity in MHz/Torr
         % out line intensity in cm-1/atm
         mb2torr = 0.750062;
         atm2torr = 1013 * mb2torr; % 1013 mb * mb2torr 
         [c,units]=constants('light');
         c = c * 1e2; %m/s -> cm/s
         out = inp * atm2torr;    % MHz/Torr -> MHz/atm
         out = (out * 1e6) / c; % MHz/atm -> cm-1*atm-1 
         units = 'cm-1/atm';
         
    otherwise
        
        disp('Unknown convertion. Sorry.');
        out = NaN;
        units = '';
        
end


return
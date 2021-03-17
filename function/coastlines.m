% Coastlines retrieves coast lines (lat,lon) for a variety of geographic
% areas and resolutions.
% Usage:
%      [nargout] = coastlines(geoarea);
% Input:
%      geoarea: a string with the desidered geographic area/resolution
% Output:
%      It depends on the geoarea
% Examples:
%      [COAST,POLBN,RIVRS] = coastlines('worldmaplo');
%      [COAST] = coastlines('gshhs_c'); 
%      [COAST] = coastlines('gshhs_l');
%      [COAST] = coastlines('gshhs_i'); 
%      [COAST] = coastlines('gshhs_h');
%      [COAST,POLBN] = coastlines('mediterraneo');
%      [COAST] = coastlines('italia');
%      [COAST] = coastlines('regioni');
%      [COAST] = coastlines('province');
%
% Nico, Jen 2008
%
% History:
% 2008/01/24 - First version

function [varargout] = coastlines(geoarea)

Matoolsfolder = '/Users/Nico/Matools/';

switch geoarea
    case 'worldmaplo'
         [COAST,POLBN,RIVRS] = worldmap; 
         % Non ho mai capito perch?, ma questa mappa sballa di 0.2 deg
         % circa (su tutto il globo)
         dlon = 0.2; 
         COAST.lon = COAST.lon - dlon; 
         POLBN.lon = POLBN.lon - dlon;
         RIVRS.lon = RIVRS.lon - dlon;
         varargout{1} = COAST;
         varargout{2} = POLBN;
         varargout{3} = RIVRS;
    case {'gshhs_c','gshhs_l','gshhs_i','gshhs_h'}
         OUT = load([Matoolsfolder 'FROMDCT/' geoarea '.mat']);
         varargout{1} = OUT;
    case 'mediterraneo'
         load Mediterraneo
         COAST.lon = Italia.lon; COAST.lat = Italia.lat; 
         POLBN.lon = [Abruzzo.lon; NaN; Lazio.lon];
         POLBN.lat = [Abruzzo.lat; NaN; Lazio.lat];
         varargout{1} = COAST;
         varargout{2} = POLBN;
    case 'italia'
         OUT = load([Matoolsfolder 'REGIONI\MFILES\Italia.mat']);
         varargout{1} = OUT;
    case 'regioni'
         OUT = load([Matoolsfolder 'REGIONI\MFILES\Regioni.mat']);
         varargout{1} = OUT;
    case 'province'
         OUT = load([Matoolsfolder 'REGIONI\MFILES\Province.mat']);
         varargout{1} = OUT;
    otherwise
         fprintf(1,'coastlines: unknown option %s.\n',geoarea);
end

return

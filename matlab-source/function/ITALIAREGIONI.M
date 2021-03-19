% This function plots boudaries for one or all 
% the italian regions.
% Usage:
%       ItaliaRegioni(Region,hf,limits);
% Input:
%       Region: string stating the wanted region
%               'All' would plot all of them
% Optinal Input:
%           hf: Figure handle for overplot
%       limits: four-component vector [LONmin LONmax LATmin LATmax]
%               default is [6 20 34 48]
% Output:
%       A Matlab Figure
%
% Nico, Sept, 2003
% Last update: 2003/24/09
%
% NB: Funzione ormai obsoleta! Vedi map_regioni e map_province. Nico 2008/01/24

function ItaliaRegioni(Region,hf,limits)


load /home/niko/Matools/REGIONI/MFILES/REGIONI.mat

if strcmp(Region,'All');

 Region = ...
 ['Abruzzo      ';...
  'Basilicata   ';... 
  'Calabria     ';... 
  'Campania     ';...  
  'EmiliaRomagna';...
  'Friuli       ';...
  'Lazio        ';...  
  'Liguria      ';...
  'Lombardia    ';... 
  'Marche       ';...
  'Molise       ';...
  'Piemonte     ';...
  'Puglia       ';...
  'Sardegna     ';...
  'Sicilia      ';...
  'Toscana      ';...
  'Trentino     ';...
  'Umbria       ';...
  'ValledAosta  ';...
  'Veneto       '];

end


if nargin>1; 
    figure(hf); 
else
    figure
end


for ir = 1:length(Region(:,1));
    
    eval(['plot(' deblank(Region(ir,:)) '.lon,' deblank(Region(ir,:)) '.lat,''k'')']);
    xlabel('Longitude');
    ylabel('Latitude');
    hold on;
        
end


if nargin>2; axis(limits); end;


return

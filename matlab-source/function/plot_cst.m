function [lon,lat] = plot_cst(lonb,latb,res,pflag)

%
% [lon,lat] = function plot_cst(lonb,latb,res,pflag)
%
% plot_cst.m returns coastline data and makes a simple line plot of 
% the data.  Use the mapping toolbox functions for a better looking 
% plot (e.g. using patchesm, displaym, etc.).
% 
% Inputs:
%  lonb,latb : Longitude and Latitude boundaries (2x1 each) of plot (deg)
%      both default to [-inf inf] if not input.
%
%  res: string referring to which resolution to use:
%     'c': crude (25km)
%     'l': low (5.0km) 
%     'i': intermediate (1.0km) or 
%     'h': high (0.2km)
%     defaults to 'l' if not input.
%  pflag = 0 for no plot,
%        = 1 for simple line plot
%        = 2 for a somwehat nicer plot with US state and country
%            borders, large rivers, etc.  requires mapping toolbox.
%        defaults to 2 if not input
%
% Output:
%  plots the data using the standard plot command (e.g. not plotm)
%  and returns lon and lat vectors (degrees)
%
% Calling Example:
%  [lon,lat] = plot_cst([-94 -81.6],[41.5 47.2],'i',2);
%
% Coastline data taken from the GSHHS: high resolution shoreline data from 
%  the National Geophysical Data Center
%  (http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html)
%  and stored in .mat files
%
% ****Bugs/Features****
% * The longitudes in the GSHHS data sets range from 0 to 360 degrees
%   (e.g. _not_ -180 to 180).  This rountine currently subtracts 360 from 
%   these values to get normal values for the US region.  A more complete
%   fix will be included in the next version.
% * There is a wrap-around line from 360 to 0 degrees longitude in 
%   Antarctica which needs to be removed.  This is done by inserting
%   a NaN between the points where diff(lon) is ~360.
%
%
% DCT 1-29-99
%

% set defaults
if exist('lonb')==0;lonb = [-inf inf];end
if exist('latb')==0;latb = [-inf inf];end
if exist('res','var')==0;res = 'l';end
if exist('pflag')==0;pflag = 2;end

% edit this line to point to the directory where the .mat files are
%data_path = '/home3/davet/mfiles/coastlines/GSHHS/';
data_path = './';

if 1 % load from .mat file
  eval(['load ' data_path 'gshhs_' res '.mat'])
else % load from ascii file
 eval(['load ' data_path 'gshhs_' res '.asc'])
 eval(['lon=gshhs_' res '(:,1);' ])
 eval(['lat=gshhs_' res '(:,2);' ])
 eval(['clear gshhs_' res])
end

% weird longitude values workaround
lon = lon-360;

% select desired region
ind = find((lon >= min(lonb)-2 & lon <= max(lonb)+2 & ...
            lat >= min(latb)-2 & lat <= max(latb))+2 | ...
            isnan(lat) == 1);
lon = lon(ind);lat = lat(ind);

% get rid of wrap-arounds
ind=find(diff(lon) > 359);
if length(ind) >= 1
  for i = 1:length(ind)
    lat = [lat(1:ind(i)) ; NaN ; lat(ind(i)+1:length(lat))];
    lon = [lon(1:ind(i)) ; NaN ; lon(ind(i)+1:length(lon))];
  end
end

[lat,lon]=maptriml(lat,lon,[min(latb) max(latb)],[min(lonb) max(lonb)]);

% plot

if pflag == 1
 
 plot(lon,lat,'Color',[.8 .8 .8],'Tag','cstplot')
 title('National Geophysical Data Center coastline data')
% axis([min(lonb) max(lonb) min(latb) max(latb)])
% if max(axis)==inf;axis auto;end
 grid on

elseif pflag == 2
 
 hold on

 % state borders
 load usahi
 for i = 1:length(stateline)
   ln=stateline(i).long;lt=stateline(i).lat;
   [lt1,ln1]=maptriml(lt,ln,[-inf inf],[-20 180]);
   [lt2,ln2]=maptriml(lt,ln,[-inf inf],[-180 -20]);
   ln=[ln1-360;ln2];lt=[lt1;lt2];
   [lt,ln]=maptriml(lt,ln,[min(latb) max(latb)],[min(lonb) max(lonb)]);
   plot(ln,lt,'Color',[.85 .5 .55],'Tag','cstplot')
 end

 % rivers, etc
 load worldlo
 ln=DNline(1).long;lt=DNline(1).lat;
 [lt1,ln1]=maptriml(lt,ln,[-inf inf],[-20 180]);
 [lt2,ln2]=maptriml(lt,ln,[-inf inf],[-180 -20]);
 ln=[ln1-360;ln2];lt=[lt1;lt2];
 [lt,ln]=maptriml(lt,ln,[min(latb) max(latb)],[min(lonb) max(lonb)]);
 plot(ln,lt,'Color',[.5 .85 .85],'Tag','cstplot')

 % political boundaries
 ln=POline(1).long;lt=POline(1).lat;
 [lt1,ln1]=maptriml(lt,ln,[-inf inf],[-20 180]);
 [lt2,ln2]=maptriml(lt,ln,[-inf inf],[-180 -20]);
 ln=[ln1-360;ln2];lt=[lt1;lt2];
 [lt,ln]=maptriml(lt,ln,[min(latb) max(latb)],[min(lonb) max(lonb)]);
 plot(ln,lt,'Color',[.85 .85 .85],'Tag','cstplot')
 
 % coastlines
 plot(lon,lat,'Color',[.6 .6 .6],'Tag','cstplot')
 title('National Geophysical Data Center coastline data')
% axis([min(lonb) max(lonb) min(latb) max(latb)])
% if max(axis)==inf;axis auto;end
 
 grid on
end  


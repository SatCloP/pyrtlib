% This function extracts a time series of a variable from a set of ARM netcdf files
%
% Usage:
%       [TSRS,OUT] = extracttimeseries(datadir,fields,wild)
% Inputs:
%        datadir: string stating the folder where cdf files are stored
%         fields: name of the needed variable(s) 
%           wild: wild card string
% Outputs:
%       TSRS: structure containing JT and all the fields asked in input by "fields"
%        OUT: structure with one component for each file
% Example:
%       [TSRS,OUT] = extracttimeseries('/home/niko/dati/RAINGB/SMOS/',['precip'; 'rh    ']);
%
% History:
%       2004/02/09 - add wild card in input
%       2007/05/26 - adapted to the new netcdf reader: nc_varget
%       2016/05/09 - adapted to the new netcdf reader: ncread

function [TSRS,OUT] = extracttimeseries(datadir,fields,wild)

% if not(exist('nc_varget')); % add ncdf tools to path
%    addpath('C:\Programmi\MATLAB71\toolbox\snctools\');
%    addpath('C:\Programmi\MATLAB71\toolbox\mexnc\');    
% end

if nargin < 3
    wild = '';
end

% 1) One structure with many components, one for each file
OUT = struct('JT',[],'YYYY',[]);
for nfld = 1:length(fields(:,1))
    eval(['OUT.' deblank(fields(nfld,:)) '=[];' ]);
end

% 2) Only one time series for all files
TSRS.JT = []; TSRS.YYYY = [];
for nfld = 1:length(fields(:,1))
    eval(['TSRS.' deblank(fields(nfld,:)) '=[];' ]);
end


listofiles = dirsort([datadir wild '*.cdf'],'name');

for nfile = 1:length(listofiles)
    
    ncfile = listofiles(nfile).name;
    fprintf(1,'Adding %s...\n',ncfile);
    
%     base_time = nc_varget([datadir ncfile],'base_time');
%     time_offset = nc_varget([datadir ncfile],'time_offset');
    base_time = ncread([datadir ncfile],'base_time');
    time_offset = ncread([datadir ncfile],'time_offset');
    base_time = double(base_time);
    [timestr,julday,datetime] = computertime(base_time + time_offset);
    
    

    % 1)
    OUT(nfile).JT = julday;    
    OUT(nfile).YYYY = datetime(:,1);
    % 2)
    sizefld = size(julday);
    if sizefld(1)==1 %row vector
        TSRS.JT = [TSRS.JT julday];
        TSRS.YYYY = [TSRS.YYYY datetime(1,:)];
    else             % colon vector
        TSRS.JT = [TSRS.JT; julday];
        TSRS.YYYY = [TSRS.YYYY; datetime(:,1)];
    end

    for nfld = 1:length(fields(:,1))
        % 1)
        fld_name = deblank(fields(nfld,:));
        %fld_data = nc_varget([datadir ncfile],fld_name);
        fld_data = ncread([datadir ncfile],fld_name);
        eval(['OUT(nfile).' fld_name ' = fld_data;' ]);
        % 2)
        sizefld = size(fld_data);
        if sizefld(1) == length(julday)
           % if nrow==ntime, thus colon vector
           eval(['TSRS.' fld_name '=[TSRS.' fld_name '; fld_data];' ]);
        else
           % if ncol==ntime, thus row vector
           eval(['TSRS.' fld_name '=[TSRS.' fld_name ' fld_data];' ]);           
        end
    end
                
end


return
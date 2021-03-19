
function [theResult] = rd_netcdf(theNetCDFFile)

%
% read in variables from a NetCDF files 
%

% get filename if not input
if nargin == 0
  [filename,pathname] = uigetfile;
  theNetCDFFile = [pathname filename];
end

% open netcdf files and get file info
f = netcdf(theNetCDFFile, 'nowrite');
if isempty(f);return;end

% get and store variable names
variables = ncnames(var(f));
theResult.variables = variables;

% loop over variable names, reading them in and assigning them to 
% the data structure if they are not the radiances.

for i = 1:length(variables)
  % if there is a '-' in the variable name, don't read it in
  if isempty(find(variables{i} == '-')) 
    eval(['theResult.' variables{i} ' = f{variables{i}}(:);'])
  end
end

close(f)

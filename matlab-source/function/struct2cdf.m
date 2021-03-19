
function rc = struct2cdf(dat,fileName);

%
% function rc = struct2cdf(dat,fileName);
%
% --------------------------------------------------------------------------
%
% NAME: struct2cdf
%
% PURPOSE: Write data in input structure "dat" to the NetCDF file "fileName".
%          Hopefully makes writing CDF files more transparent to the user.
%
% COPYRIGHT INFORMATION:  
%	  (C) COPYRIGHT UW/SSEC ALL RIGHTS RESERVED, 2001
%	      Space Science & Engineering Center
%	      University of Wisconsin - Madison  
% INPUTS:
%   dat    is a structure containing the data and attributes to be written to the 
%            NetCDF file.  dat should contain the following fields:
%              .dim         contains all of the dimensions (Required)
%              .globalatt   contains any global attributes (Optional)
%            and for each variable, <variable_name>, 
%              .<variable_name>.dim   is a cell array defining the dimensions 
%                                    of this variable (Required).  These strings must 
%                                    match field names of .dim.
%              .<variable_name>.type   defines the data type, using strings defined 
%                                     with "nctype". (Optional - defaults to "ncfloat").
%              .<variable_name>.val   contains the actual data. (Required)
%              .<variable_name>.att   is a structure containing any attributes for this 
%                                  variable (Optional).
%
%   fileName  is the name of the NetCDF file to be created.
%
% OUTPUTS:
%      rc       return code: ==1 if the file was created sucessfully, ==0, if not
%
% OUTPUT FILES: The file <fileName> is created in the user's local directory.
%                If this file already exists, it may be overwritten.
%
% EXAMPLE USAGE:
%   % create input structure
%   clear dat
%   dat.dim.time = 24;
%   dat.dim.other = 5;
%   dat.hourOfDay.dim = {'time'};
%   dat.hourOfDay.type = {'ncint'};
%   dat.hourOfDay.val = (0:23)';
%   dat.hourOfDay.att.longname = 'Hour of the Day';
%   dat.hourOfDay.att.units = 'hour';
%   dat.dummyData.dim = {'time','other'};
%   dat.dummyData.val = zeros(24,5);
%   dat.globalatt.author = {'me'};
%   dat.globalatt.created_on = datestr(now);
%   % write the file
%   rc = struct2cdf(dat,'sample_file.cdf');
%
% INTERNAL or EXTERNAL CALLS:
%   Uses the USGS Matlab/NetCDF toolbox for writing the NetCDF file.
%   (http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html)
%
% RESTRICTIONS: ??
%
% RELATED: netcdf, rd_netcdf
%
% MODIFICATION HISTORY: 
%	  Written by:	 Dave Tobin,  02-05-2001 (UW-SSEC)
%
% REVISION INFORMATION (RCS Keyword):
%   $Id: struct2cdf.m,v 1.2 2001/05/02 12:28:02 davet Exp $
%
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%                REVISION INFORMATION (RCS Keyword variable)
%------------------------------------------------------------------------------
rcs_id =  '$Id: struct2cdf.m,v 1.2 2001/05/02 12:28:02 davet Exp $' ;

% Initialize warning mode
warning backtrace

% Initialize return code to valid (1)
rc = 1;

% Set number of required input variables (Set manually)
nargin_req = 2;

% Test for a valid number of input arguments
if nargin < nargin_req
  warn_msg = ['Not enough input arguments to ' mfilename ' ... aborting'];
  rc = 0;
elseif ~isfield(dat,'dim');
  warn_msg = ['dimensions not defined in dat.dim ... aborting'];
  rc = 0;
end

% Display warnings and abort if necessary
if rc == 0
  eval(['help ' mfilename])
  warning(warn_msg)
  return
end


try

% --------------------------------------
% Open the output file
% --------------------------------------
ncquiet;
nc = netcdf(fileName,'clobber');

% -----------------------
% Define dimensions
% -----------------------
dimensions = fieldnames(dat.dim);
ndimensions = length(dimensions);
for i = 1:ndimensions
  nc(char(dimensions(i))) = getfield(dat.dim,char(dimensions(i)));
end

% -------------------------------------------------------------------
% Write Global Attributes, if present in input structure
% -------------------------------------------------------------------
if isfield(dat,'globalatt')
  globalAttNames = fieldnames(dat.globalatt);
  nglobalAtts = length(globalAttNames);
  for i = 1:nglobalAtts
    eval(['nc.' char(globalAttNames(i)) ' = ''' ...
		char(getfield(dat.globalatt,char(globalAttNames(i)))) ''';']);
  end
else
  nglobalAtts = 0;
end

% -------------------------
% Loop over variables
% -------------------------
variables = fieldnames(dat);
nvariables = length(variables);
for i = 1:nvariables

  this_variable_name = char(variables(i));
  this_variable = getfield(dat,this_variable_name);

  % Skip variables which do not have the dimension field "dim".
  if isfield(this_variable,'dim')

    % Define dimensions
    if ~isfield(this_variable,'type')
      this_variable.type = {'ncfloat'};
    end
    tmp = [];
    for j = 1:length(this_variable.dim)
    tmp = [tmp '''' char(this_variable.dim(j)) '''' ','];
    end
    tmp = tmp(1:end-1);
    eval(['nc{''' this_variable_name '''} = ' char(this_variable.type) '({' tmp '});']);

    % Write data
    nc{this_variable_name}(:) = this_variable.val;

    % Write attributes for this variable if present
    if isfield(this_variable,'att')
      attNames = fieldnames(this_variable.att);
      nAtts = length(attNames);
      for i = 1:nAtts
        eval(['nc{''' this_variable_name '''}.' char(attNames(i)) ' = ''' ...
		char(getfield(this_variable.att,char(attNames(i)))) ''';']);
      end
    end
  end
end

%--------------------
% Close the file
%--------------------
close(nc);


catch

  rc = 0;
  warning(['An error occured when creating : ' fileName ' ... aborting']);
  return

end

%------------------------------------------------------------------------------
%
% EXPANDED REVISION INFORMATION (RCS Keywords):
%
%   $Id: struct2cdf.m,v 1.2 2001/05/02 12:28:02 davet Exp $
%   $Author: davet $
%   $Log: struct2cdf.m,v $
%   Revision 1.2  2001/05/02 12:28:02  davet
%   *** empty log message ***
%
%   Revision 1.1  2001/05/02 12:26:03  davet
%   *** empty log message ***
%
%   $Locker:  $ 
%
%------------------------------------------------------------------------------

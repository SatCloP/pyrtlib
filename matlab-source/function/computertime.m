% This function computes the date that corrispondes
% to input (scalar or vector) seconds from 1/1/1970.
% If you want a different "0" time instead of 1/1/1970,
% have to give in input the date of your own "0" time.
%
% Es:
%    [timestr,julday,datetime]=computertime(secs)
%    or
%    [timestr,julday,datetime]=computertime(secs,date0)
%
% Inputs:
%    secs: Seconds from 00:00:00 of date0 (Default: 1/1/1970) 
%   date0: Date of your own "0" time. 
%          Three component vector [YYYY MM dd] (Es: [1987 1 1])
%
% Outputs:
%   timestr: A string with complete date and time
%    julday: Decimal julian day
%  datetime: Vector (or Matrix) containing YY,MM,dd,hh,mm,ss
%
% History
%   2007/12/18 - Fixed a bug to force dimension consistency

function [daystr,julday,datetime] = computertime(secs,date0)

  if nargin < 2
    day = secs/(86400) + datenum(1970,1,1);
  else
    day = secs/(86400) + datenum(date0(1),date0(2),date0(3));
  end
  
  daystr = datestr(day);
    
  thatyear = str2num(datestr(day,10));
  
  day = reshape(day,size(thatyear));
  
  julday = day-datenum(thatyear,1,1)+1;

  datetime = datevec(day);
  
return
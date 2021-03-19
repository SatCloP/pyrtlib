% FUNCTION THAT CALCULATES DATE (MONTH, DAY AND ALSO HOURS,  
% MINUTES AND SECONDS) GIVING IN INPUT YEAR AND JULIAN DAY. 
% IT ACTUALLY TAKES CARE OF LEAP YEARS.
% NB: IT GIVES BACK THE UTC DATE (NOT THE LOCAL ONE)
%
%  ES:
%      daystr=jul2date(YYYY,julday);
%      [MM,dd,hh,mm,ss]=jul2date(YYYY,julday);
%
% Nico
% LAST UPDATE: 4 MAR 2012 - fixed a bug in the consistency check (julday > sum(dforM)+1)


function [MM,dd,hh,mm,ss]=jul2date(YYYY,julday);

% CHECKING INPUT
dforM=[31 28 31 30 31 30 31 31 30 31 30 31];
if (mod(YYYY,4)==0)
  dforM(2)=29;   
end
if any( julday > sum(dforM)+1 )
   fprintf(2,'\n Are you sure of this julian day: %f ',julday)
   fprintf(2,'\n Mmmh... Is greater then the whole year.\n')
   return
end

if nargout<=1
   MM=datestr(datenum(YYYY,1,1) + julday - 1,0);
else   
   [DUM,MM,dd,hh,mm,ss] = datevec(datenum(YYYY,1,1) + julday - 1);
end

return

% FUNCTION THAT GIVES THE DECIMAL JULIAN DAY, GIVING IN INPUT THE DATE 
% (YEAR, MONTH, DAY AND OPTIONALLY HOURS, MINUTES AND SECONDS).
% IT TAKES CARE OF LEAP YEARS.
% NB: IT WANTS THE UTC DATE (NOT THE LOCAL ONE)
%
% Es: [julday]=date2jul(YYYY,MM,dd,hh,mm,ss);
%
% Nico
% LAST UPDATE: 9 MAR 2001

function [julday]=date2jul(YY,MM,dd,hh,mm,ss);

switch nargin
case 3
   hh=0; mm=0; ss=0;
case 4
   mm=0; ss=0;   
case 5
   ss=0;   
end

% CHECKING INPUT %%%%%%%
%dforM=[31 28 31 30 31 30 31 31 30 31 30 31];
%if (mod(YY,4)==0)
%  dforM(2)=29;   
%end
%if MM>12 | dd>dforM(MM) | hh>24 | mm>60 | ss>60
%   fprintf(2,'\n Are you sure of this date? YY:%d MM:%d dd:%d hh:%d mm:%d ss:%d',YY,MM,dd,hh,mm,ss)
%   fprintf(2,'\n Mmmh... Something seems to be wrong.\n')
%   return
%end
%%%%%%%%%%%%%%%%%%%%%%%%

julday = datenum(YY,MM,dd,hh,mm,ss) - datenum(YY,1,1) + 1 ;

return
   
   
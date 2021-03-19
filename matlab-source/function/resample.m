% Function that resample a vector y defined on x,
% on a new xn, based on the nearest values within
% a range given in input wd.
% If no points is found within wd the resampled 
% function will have a "missing value" flag,
% also specified in input, mv.
%
% Es:
%    yn=resample(x,y,xn,wd,mv);
% Inputs:
%        x  Initial grid
%        y  Initial vector
%        xn Resempling grid
%        wd Range Width (units as x)
%        mv Missing value flag
%
% Nico, 10/2000

function yn=resample(x,y,xn,wd,mv);

 if length(x)~=length(y); 
    fprintf(1,'Different size for x and y'); 
    yn=[]; 
    return; 
 end;

 lxn=length(xn);
 
 for i=1:lxn 
    
   [dif,nearst] = min( abs(x-xn(i)) );
   %fprintf(1,'i = %d \n',i); 
 
   if dif < wd
      yn(i) = y(nearst);
   else
      yn(i) = mv;    
   end

 end
 
return
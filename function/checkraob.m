% This function checks the input raob measures (z,p,tk,rh)
% to see if they are reasonable or not.
% It gives back the same quantities, but just for those 
% level that passed the quality control.
%
% Quality control applied:
% 1) Check for Relative Umidity out of range ([0 1])
% 2) Check for crazy Temperature ([0 340]K)
% 3) Check for crazy Height ([0 50]Km)
% 4) Check for crazy Pressure ([0 1500]mb)
% 5) Check for monothonic Pressure
% 6) Check for Temperature Gradient (10K)
%
% Es:
%    [z,p,tk,rh]=checkraob(z,p,tk,rh);
%
% Inputs:
%        z   Height (Km)
%        p   Pressure (mb) 
%        tk  Temperature (K)
%        rh  Relative Umidity (%/100)
%
% Nico, 10/2000
% History:
%              04 Mar 2004 - Add NaN removal
% Last update: 24 Mar 2004 - Debuged problem with line vector

function [z,p,tk,rh]=checkraob(z,p,tk,rh)

% 1) Check for Relative Umidity out of range
badrh=find(rh > 1 | rh < 0);

% 2) Check for crazy Temperature
badtk1=find(tk > 340 | tk < 0);

% 3) Check for crazy height 
badz=find(z > 50 | z < 0);

% 4) Check for crazy Pressure
badp1=find(p > 1500 | p < 0);

% 5) Check for monothonic pressure
dp=diff(p);
badp2=find(dp > 0);
badp2=badp2+1;

% 6) Check for Temperature Gradient
dtk=diff(tk);
badtk2=find(abs(dtk) > 10);

% Purging bad levels
%bad=union([badrh; badtk1; badtk2],[badz; badp1; badp2]);
bad=union(union(badtk2,union(badrh,badtk1)),union(badz,union(badp1,badp2)));
z(bad)=[];
p(bad)=[];
tk(bad)=[];
rh(bad)=[];

% 7) Check for NaNs
nl = size(z);
if nl(1)==1 % row vector
   [inan,jnan] = find(isnan([z' p' tk' rh']));
else        % colon vector
   [inan,jnan] = find(isnan([z p tk rh]));
end    
z(inan)=[];
p(inan)=[];
tk(inan)=[];
rh(inan)=[];

return
% Function that gives back the local time if provided
% UTC and time shift (+/-).
% I know it's obvious but...still you can save 3 lines!
%
% Es:
%    [loctime]=localtime(utc,shift)
%    or
%    [utc]=localtime(loctime,shift)  
%
% Nico,2000

function [loctime]=localtime(utc,shift)

loctime=utc+shift;
if loctime<0
   loctime=loctime+24;
elseif loctime>24
   loctime=loctime-24;
end

return
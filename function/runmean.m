% This function computes a running mean of the input vector.
% A running mean doesn't reduce the sample size; if you want 
% to reduce your sample size check out function boxmean.
%
% Es: 
%        [AVR,SGM,tav]=runmean(V,t,hm)
% Inputs:
%        V: vector to average
%        t: time (same size of V)
%       hm: time interval to divide the time range
% N.B.: It doesn't matter what units are used for time
%       but have to be the same either for time and hm
%
% Outputs:
%      AVR: average of V inside each sigle bin of width hm
%      SGM: standard deviation in each bin
%
% Note: It doesn't reduce the sample size.
%       If it's needed a component running mean 
%       (instead of time), it can be done using 
%       [AVR,SGM,tav]=runmean(V,1:length(V),nc);
%       where nc is the number of components to average.
%
% Nico, 2001 (OLDER VERSION:2000)

function [AVR,SGM]=runmean(V,t,hm)

% Make sure time is increasing
[t,indx] = sort(t);
V = V(indx);
clear indx

for i = 1:length(V)
   
   strt = t(i)-hm/2;
   stop = t(i)+hm/2;
   
   indx = find (t>strt & t<stop);
   
   if isempty(indx)
     AVR(i) = NaN;
     SGM(i) = NaN;
   else      
     AVR(i) = mean(V(indx));
     SGM(i) = std (V(indx));
   end
   clear indx

end

% Purging NaN
indx = find(isnan(AVR));
AVR(indx) = [];
SGM(indx) = [];

return
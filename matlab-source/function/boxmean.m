% This function computes a box mean of the input vector.
% A box mean reduces the sample size; if you don't want 
% this check out function runmean.
%
% Es: 
%        [AVR,SGM,tav]=boxmean(V,t,hm)
% Inputs:
%        V: vector to average
%        t: time (same size of V)
%       hm: time interval to divide the time range
% Outputs:
%      AVR: average of V inside each sigle bin of width hm
%      SGM: standard deviation in each bin
%      tav: time average inside each single bin
%
% Note: It reduces the sample size.
%       If it's needed a component box mean 
%       (instead of time), it can be done using 
%       [AVR,SGM,tav]=boxmean(V,1:length(V),nc);
%       where nc is the number of components to average.
%
% N.B.: It doesn't matter what units are used for time
%       but have to be the same either for time and hm
%
% Nico, 2001 (OLDER VERSION:1999)

function [AVR,SGM,tav]=boxmean(V,t,hm)

% Make sure time is increasing
[t,indx] = sort(t);
V = V(indx);
clear indx

% Create time bins
nbins = floor( (t(end)-t(1)) / hm );

% Filling the bins
for ibin = 1:nbins
   
   strt = t(1)+(ibin-1)*hm;
   stop = t(1)+(ibin)*hm;
   
   indx = find (t>strt & t<stop);
   
   if isempty(indx)
     AVR(ibin) = NaN;
     SGM(ibin) = NaN;
     tav(ibin) = NaN;
   else      
     AVR(ibin) = mean(V(indx));
     SGM(ibin) = std (V(indx));
     tav(ibin) = mean(t(indx)); 
   end
   clear indx

end

% Purging NaN
indx = find(isnan(tav));
tav(indx) = [];
AVR(indx) = [];
SGM(indx) = [];

return
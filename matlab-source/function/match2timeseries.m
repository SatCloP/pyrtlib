% Usage:
% [var2attime1,missvalperc] = match2timeseries(time1,time2,var2,int,missval)
%
% 2010/08/23 - Added the missing value option (default -999.99)

function [var2attime1,missvalperc] = match2timeseries(time1,time2,var2,int,missval)

if nargin < 5
   missval = -999.99;
end

nt = length(time1);
var2attime1 = NaN * ones(size(time1));
missvalperc = NaN * ones(size(time1));

for it = 1:nt
    
    indx = find( abs( time2-time1(it) ) <= int );
    if ~isempty(indx)
%        var2attime1(it) = nanmean( var2(indx) );
%        missvalperc(it) = length( find(var2(indx)==missval) ) / length(indx); % ranges from 0 (none missing) to 1 (all missing)
        dum = var2(indx); 
        indxmiss = find( dum == missval ); dum( indxmiss ) = NaN;
        var2attime1(it) = nanmean( dum );
        missvalperc(it) = length( indxmiss ) / length(indx); % ranges from 0 (none missing) to 1 (all missing)
    end

end


return

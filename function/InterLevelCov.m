% This function computes the inter-level error covariance that can be used
% as an estimation of vertical resolution, as described in: 
% Liljegren, J.C., S. A. Boukabara, K. Cady-Pereira, and S.A. Clough, The Effect of the 
%   Half-Width of the 22-GHz Water Vapor Line on Retrievals of Temperature and Water Vapor 
%   Profiles with a Twelve-Channel Microwave Radiometer, Trans. Geosci. Remote Sensing, 2005.
% Gueldner, J., and D., Spaenkuch, Remote sensing of the thermodynamic state of the atmospheric
%   boundary layer by microwave radiometry, J. Atmos and Ocean. Tech., vol. 18, pp. 925-933, 2001.
% Smith, W.L., W.F. Feltz, R.O. Knuteson, H.R. Revercomb, H.B. Howell, and H.H. Wolf, 
%   The retrieval of planetary boundary layer structure using ground-based infrared spectral 
%   radiance measurements, J. Atmos. Ocean. Tech., vol. 16, pp. 323-333, 1999.
%
% Usage:
%       [R,C] = InterLevelCov(Z,YE,YR,removebias);
% Input:
%       Z: heigth (nlev) [whatever units are given for Z, R would have the same]
%      YE: estimated profiles (nlev x nprof) [same units as for YR]
%      YR: radiosonde profiles (nlev x nprof) [same units as for YE]
% removebias: if 1 removes bias between YE and YR (which may affect the resolution) 
% Output:
%       R: resolution (nlev) [same units as for Z]
%       C: inter-level error covariance (nlev x nlev) [dimensionless]
% Example:
%       [R,C] = InterLevelCov(Z,YE,YR);
%
% Nico, Aug 2004
%
% History:
% 2004/08/04 - First version
% 2005/06/14 - Add help
% 2005/06/22 - Fast algorithm
% 2008/03/09 - Added help for removebias

function [R,C] = InterLevelCov(Z,YE,YR,removebias)

nlev = length(Z);
nprof = length(YE(1,:));
C = zeros(nlev,nlev);
R = zeros(nlev,1);
thr = 0.6;

if removebias
   bias = mean(YE-YR,2); 
   YE = YE - bias*ones(1,nprof);
end

NUM = ((YE-YR)*(YE-YR)');
DEN = sqrt( sum(((YE-YR).^2)')' * sum(((YE-YR).^2)') );
C = NUM ./ DEN;

for nl1 = 1:nlev
    
     %for nl2 = 1:nlev
     %     NOM = sum( (YE(nl1,:)-YR(nl1,:)) .* (YE(nl2,:)-YR(nl2,:)) );
     %     DEN = sqrt( sum( (YE(nl1,:)-YR(nl1,:)).^2 ) * sum( (YE(nl2,:)-YR(nl2,:)).^2 ) );
     %     C(nl1,nl2) = NOM./DEN;
     %end
     
     indx = find(C(nl1,:)>thr);
     
     if isempty(indx)
         R(nl1) = NaN;
     else
        jumps = find(diff(indx)>1);
        if ~isempty(jumps);
           jumps = [0 jumps length(indx)];
           for ig = 1:length(jumps)-1; groups(ig,1)=jumps(ig)+1; groups(ig,2)=jumps(ig+1); end; 
           group = find( nl1>=indx(groups(:,1)) & nl1<=indx(groups(:,2)) );
           indx = indx(groups(group,1)):indx(groups(group,2));
           clear groups group jumps
        end
        R(nl1) = Z(indx(end)) - Z(indx(1)); 
     end

 end


return
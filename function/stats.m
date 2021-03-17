% stats.m compute with one command median, std, max, min and mean of input vector
%
% Inputs:
%        V: N-dimension Vector 
% Optional Input:
%     indx: index of V components on which reduce the statistical analysis
% Output:
%     MEDN: Median of V
%     STDV: Standard deviation of V
%     MAXV: Maximum value of V
%     MINV: Minimum value of V
%     MEAN: Mean value of V
%     RMSE: rmse value of V
%
% Es:
%    [MEDN,STDV,MAXV,MINV,MEAN,RMSE]=stats(V)
%    [MEDN,STDV,MAXV,MINV,MEAN,RMSE]=stats(V,indx)
%
% 2008/06/05 - Add RMSE

function [MEDN,STDV,MAXV,MINV,MEAN,RMSE]=stats(V,indx)

if nargin < 2
  MEDN = median(V);
  STDV = std(V);
  MAXV = max(V);
  MINV = min(V);   
  MEAN = mean(V);   
  RMSE = rms(V);   
else
  MEDN = median(V(indx));
  STDV = std(V(indx));
  MAXV = max(V(indx));
  MINV = min(V(indx));
  MEAN = mean(V(indx));   
  RMSE = rms(V(indx));   
end

return
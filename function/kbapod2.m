
function apod = kbapod2(N, iMD, k)

%
% function apod = kbapod2(N, iMD, k)
%
% compute Kaiser-Bessel apodization function
% inputs:
%   N:  length of apodization function (double sided ifg)
%   iMD: index to point where opd = MOPD (for single sided ifg)
%   k - Kaiser-Bessel parameter (int > 0).  DEFAULT=6
%
% output
%   apod - apodization function
%
% This is a modified version of Howard Motteler's original kbapod.m 
% In this version, the resulting function is arranged for MATLAB style
% ffts with ZPD at apod(1) and MOPD at apod(N/2) (in the same manner as 
% beer.m).
%
% DCT 2-12-98
%
if ~exist('k'); k=6; end
apod = zeros(N,1);
d = linspace(0,N/2,N/2)';
x = k * sqrt(1 - (d/iMD).^2) ;
% I(x)
r = ones(N/2,1); f = 1;
for j = 1:8;
  f = f * j;
  r = r + (x/2).^(2*j) / f^2 ;
end

% I(k)
s = 1; f = 1;
for j = 1:8;
  f = f * j;
  s = s + (k/2).^(2*j) / f^2 ;
end
c = (abs(d) <= iMD) .* r / s;
apod(1:N/2)=c;
apod(N/2+1:N) = flipud(c);


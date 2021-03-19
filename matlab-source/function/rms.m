% RMS CALCULATES THE ROOT MEAN SQUARE WITHOUT 
% REMOVING THE MEAN (AS INSTEAD STD DOES).
% IF FLAG=0 IT NORMALIZES BY N-1, OTHERWISE BY N.
%
% Es:
%     rms(x);
%     rms(x,flag);
%     rms(x,flag,dim);
%
% TYPE HELP STD FOR MORE REFERENCE.
%
% Nico 6/2000

function y = rms(x,flag,dim)

if length(x)==1; y=0; return; end;
if nargin<2, flag = 0; end
if nargin<3, 
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

if flag
  y = sqrt( sum(x.*x,dim) / size(x,dim) );
else
  y = sqrt( sum(x.*x,dim) / (size(x,dim)-1) );
end


return
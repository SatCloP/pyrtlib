% This function computes hanning apodization function.
%
% Usage:
%       apod = hanning(N,MD)
% Inputs:
%        N: length of apodization function (double sided ifg)
%       MD: index to point where opd = MOPD (for single sided ifg)
% Outputs:
%     apod: apodization function
%
% Comments: resulting function is arranged for MATLAB style
%           ffts with ZPD at apod(1) and MOPD at apod(N/2)
%

function apod=hanning(N,MD);

hanning = zeros(N,1);

for i=1:N/2
   if i <= MD
      hanning(i) = 0.5 * (1. + cos(pi*i/MD) );
   else
		hanning(i)=0;
	end
end

% The following might be the frequency correction that Ed and Joe were talking about. 
% We won't use that.
if 0  
  v0 = 1000.;
  vlaser = 15800.;
  xmax = N/vlaser;
  b = 0.030;
  z = 0.5*pi*v0*xmax*b^2;
  for i=1:N/2
     hanning(i) = hanning(i)*( 1 - (z*i/(N/2))^2/6.);
  end   
end

hanning(N/2+1:N) = flipud(hanning(1:N/2));
apod = hanning;

return


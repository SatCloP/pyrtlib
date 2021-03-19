% This function computes Beer apodization function.
%
% Usage:
%       apod = beer(N,MD)
% Inputs:
%        N: length of apodization function (double sided ifg)
%       MD: index to point where opd = MOPD (for single sided ifg)
% Outputs:
%     apod: apodization function
%
% Comments: resulting function is arranged for MATLAB style
%           ffts with ZPD at apod(1) and MOPD at apod(N/2)
%

function apod=beer(N,MD);

beer=zeros(N,1);
beer(1)=1;

for i=2:N/2
	if i <= MD
		beer(i)=(1-((i-1)/MD)^2)^2;
	else
		beer(i)=0;
	end
end

beer(N/2+1:N)=flipud(beer(1:N/2));
apod=beer;

return

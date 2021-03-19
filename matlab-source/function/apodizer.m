
function [rout] = apodizer(rin,aflag);

%
% function [rout] = apodizer(rin,aflag);
%
% apodize input sprectrunm, rin
%
% aflag :determines which apodizationfunctionto use:
% = 1, Beer
% = 2, Kaiser Bessel,k=6
% = 3, Hanning
%
% Assumes the input spectrum, rin, is minimally sampled in
% spectral space.  e.g. for an input spectrum of N points, 
% the N/2 point of the ifg is the MaxOPD point.
%
% DCT 2-18-99
%
% History
% 2005/03/16 - Added hanning apodization function (Nico)


N = length(rin);
iMD = length(rin)/2;

if aflag == 1
 apod = beer(N,iMD);
elseif aflag == 2
 apod = kbapod2(N,iMD,6);
elseif aflag == 3
 apod = hanning(N,iMD);    
end

rout = fft(ifft(rin).*apod);

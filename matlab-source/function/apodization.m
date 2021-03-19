% This function computes apodization on atmospheric emission spectra measurements.
%
% Usage:
%       [Rout] = apodization(Rin,aflag);
%
% Inputs: 
%        Rin: Input Radiance
%      aflag: determines which apodization function to use ('beer',hanning','KaiserBessel')
% Output: 
%       Rout: Apodized Radiance
%
% Assumes the input spectrum, rin, is minimally sampled in
% spectral space.  e.g. for an input spectrum of N points, 
% the N/2 point of the ifg is the MaxOPD point.
%
% DCT 2-18-99
%

function [Rout] = apodization(Rin,aflag);

N = length(Rin);
itwasodd=0;
if floor(N/2)~=N/2;
   itwasodd=1;
   Rin(N+1)=0;
   N=N+1;
end   
iMD = length(Rin)/2;

switch aflag
case 'beer'
   apod = beer(N,iMD);
case 'hanning'
   apod = hanning(N,iMD);   
case 'KaiserBessel'   
   apod = kbapod2(N,iMD,6);
end

Rout = fft(ifft(Rin).*apod);

if itwasodd; Rout(N) = []; end;

return
% This function downscale the resolution of a matrix 
% keeping the original dimension.
%
% Usage: 
%       [AMTX,STD] = averagegrid(MTX,nc,nr);
% Inputs:
%       MTX: matrix (NxM)
%        nc: colon steps to include in the average (i-nc:1:i+nc)
%        nr: row steps to include in the average (i-nr:1:i+nr)
% Outputs:
%      AMTX: downscaled matrix (NxM)
%      STDX: matrix of std over the ncxnr subgrid (NxM)
%
% Nico, Nov 2003

function [AMTX,STDX] = averagegrid(MTX,nc,nr);

[nrows,ncols] = size(MTX);
AMTX = zeros(size(MTX));
STDX = zeros(size(MTX));

for irow = 1:nrows
    for icol = 1:ncols        
        indxrow = max(1,irow-nr) : 1 : min(nrows,irow+nr);
        indxcol = max(1,icol-nc) : 1 : min(ncols,icol+nc);
        npoints = length(indxrow)*length(indxcol);
        VECT(1:npoints) = MTX(indxrow,indxcol);
        AMTX(irow,icol) = mean( VECT );
        STDX(irow,icol) = std( VECT );
        clear VECT npoints indxrow indxcol
    end
end


return
% This function gives back the nearest point of a grid to a given location 
% once the coordinates of the location and of the grid are provided.
% 
% Usage:
%    [COORD,INDX] = findnearest(lat1,lon1,LAT,LON)
% Input:
%    lat1,lon1: coordinates of a location (scalars)
%    LAT, LON : coordinates of a grid (matrices)
% Output:
%    COORD: coordinates (lat,lon) of the nearest point
%    INDX : grid indx of the nearest point (irow, jcol)
% Es:
%
% History
% 200? - First version
% 2008/01/24 - Changed algorithm

function [COORD,INDX] = findnearest(lat1,lon1,LAT,LON)
 
  %nlat = length(LAT(:,1));
  %nlon = length(LON(1,:));
  
  %LAT = LAT(1:end);
  %LON = LON(1:end);

  DKM = mydistance(lat1,lon1,LAT,LON);
    
  %[DUM,indx] = min(DKM);

  %ilon = floor(indx/nlat) + 1;
  %ilat = indx - (ilon-1)*nlat;
  %COORD = [LAT(indx) LON(indx)];
  %INDX = [ilat ilon];

  [X,i] = min(DKM);
  [Y,j] = min(X); i = i(j);
  COORD = [LAT(i,j) LON(i,j)];
  INDX = [i j];
    
return
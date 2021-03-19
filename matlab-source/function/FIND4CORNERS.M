% Trova i quattro punti piu' vicini di una griglia (regolare).
%
% Ancora da testare e migliorare....


function [lnb,ltb]=find4corners(ilon,ilat,lon,lat)

lonrng = [min(lon) max(lon)];
if ilon < lonrng(1) | ilon > lonrng(2)
   fprintf(1,'find4corners: Longitude %f outside input range [%f %f]',ilon,lonrng(1),lonrng(2));
   return
end

latrng = [min(lat) max(lat)];
if ilat < latrng(1) | ilat > latrng(2)
   fprintf(1,'find4corners: Latitude %f outside input range [%f %f]',ilat,latrng(1),latrng(2));
   return
end

[dum] = find(lon <= ilon);
lnb= [dum(end) dum(end)+1];
[dum] = find(lat <= ilat);
ltb= [dum(end) dum(end)+1];

return

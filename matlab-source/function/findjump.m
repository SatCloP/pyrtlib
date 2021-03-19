% QUESTO PROGRAMMINO MI E' SERVITO UNA VOLTA PER 
% PORTARE DEI FILE ASCII LAT-LON (ABRUZZO, LAZIO, ITALIA)
% IN FILES MATLAB

fid=fopen('contour.tmp','w');

lon=C(:,1);
lat=C(:,2);
flag=C(:,3);

for i=1:length(lon)
   
   if flag(i)==0; 
      fprintf(fid,'   NaN        NaN \n');
   end
   
   fprintf(fid,'%10.4f %10.4f \n',lon(i),lat(i));
   
end

fclose(fid);
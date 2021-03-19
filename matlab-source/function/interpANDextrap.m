% interpANDextrap allows to interpolate data points
% in a range bigger then the original, avoiding any output NaN.
% It use interp1 for linear interpolation, while extrapolate 
% linearly if necessary for points outside the original range.
% 
% Es:
%    YI=interpANDextrap(X,Y,XI,howmany)
%
% Inputs:
%    X,Y,XI The same as for interp1
%    howmany: how many points to use in the extrapolation ('all' or a number > 1)
%
% Nico, 10/18/2000

function YI=interpANDextrap(X,Y,XI,howmany)


YI=interp1(X,Y,XI);

out=isnan(YI);

indx=find(out);

NY=extrap(X,Y,XI,howmany);

YI(indx)=NY;
      
if ~isempty(find(isnan(YI)))
   fprintf(1,'interpANDextrap: Unrecoverable NaN. Sorry....\n')
end

return


function NY=extrap(X,Y,NX,hm)

[X,I]=sort(X);
Y=Y(I);

if strcmp(class(hm),'double')
   XH=X(end-hm:end);
   YH=Y(end-hm:end);
   XL=X(1:hm);
   YL=Y(1:hm);
else   
   XH=X;
   YH=Y;
   XL=X;
   YL=Y;
end

pH=polyfit(XH,YH,1);
pL=polyfit(XL,YL,1);

IH=find(NX > max(X));
IL=find(NX < min(X));

NYH=polyval(pH,NX(IH));
NYL=polyval(pL,NX(IL));

NY=[NYL NYH];

return
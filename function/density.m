% FUNCTION DENSITY
% IT COMPUTES THE 2D HISTOGRAM OF THE PREDICTED AREA
% IT 

function [npoint,Vx,Vy]=density(X,Y,stepX,stepY);

Xlim=[floor(min(X)) ceil(max(X))]
Ylim=[floor(min(Y)) ceil(max(Y))]

Vx=Xlim(1):stepX:Xlim(2);
Vy=Ylim(1):stepY:Ylim(2);

npoint=zeros(length(Vx),length(Vy));

for i=1:length(X)
   
   %fprintf('raob # %d\n',i);
   cp=[X(i) Y(i)];
   
   for h=1:length(Vx)-1
      
      for j=1:length(Vy)-1
      
         if cp(1)>=Vx(h) & cp(1)<Vx(h+1) & cp(2)>=Vy(j) & cp(2)<Vy(j+1) 
         
            npoint(h,j)=npoint(h,j)+1;
         
         end
      
      end
      
   end
   
end

% SI PUO BUTTARE

function Tauhist

nraob=[2353 1492 1179];
%nraob=[3 4 2];
ntest=12;
indx=[2 7 13 14 15 16 17 18 19 21 22 23];
Chnl=['20.375      ';'20.6        ';'21.485      ';'22.235      ';'22.985      ';'23.735      ';'23.8        ';'31.4        ';'31.65       ';'36.5        ';'89.         ';'150.        ';'183.31+/-0.5';'183.31+/-1. ';'183.31+/-3. ';'183.31+/-5. ';'183.31+/-7. ';'183.31+/-12.';'183.31+/-15.';'220.        ';'325.+/-1.   ';'325.+/-3.   ';'325.+/-8.   ';'340.        '];

nb=input('Which sample of data? (1,2,3):');


   filein=['../../data/TbTmr906_' num2str(nb) '.dat'];
   fid1=fopen(filein,'r');

   for i=1:nraob(nb)
      
      surf=carica(fid1,4,1,'%f');
      pwv(i)=surf(1);
      Psur(i,1)=surf(2);
      Tsur(i,1)=surf(3);
      RHsur(i,1)=surf(4);
      %pwv(i)=carica(fid1,1,1,'%f');
      MTb=carica(fid1,34,3,'%f');
      [Tb(i,:),Tmr(i,:),Tau(i,:)]=sel(MTb);
      clear MTb;
      clear surf;
      
   end

   for j=1:ntest
       subplot(3,4,j)
       set(gca,'Fontsize',9);
       hist(Tau(:,indx(j)))
       title(['f=' Chnl(indx(j),:) '(Np)'])
   end

clear all
fclose('all')

% FUNCTION CARICA
function MX=carica(fid,nc,nr,format);

   MX=(fscanf(fid,format,[nc,nr]))';

return

% FUNCTION SEL
function [Tb,Tmr,Tau]=sel(MTbef);

       MTb=MTbef(:,1:12);
       MTb(:,13)=(MTbef(:,19)+MTbef(:,20))/2;
       MTb(:,14)=(MTbef(:,18)+MTbef(:,21))/2;
       MTb(:,15)=(MTbef(:,17)+MTbef(:,22))/2;
       MTb(:,16)=(MTbef(:,16)+MTbef(:,23))/2;
       MTb(:,17)=(MTbef(:,15)+MTbef(:,24))/2;
       MTb(:,18)=(MTbef(:,14)+MTbef(:,25))/2;
       MTb(:,19)=(MTbef(:,13)+MTbef(:,26))/2;
       MTb(:,20)=MTbef(:,27);
       MTb(:,21)=(MTbef(:,30)+MTbef(:,31))/2;       
       MTb(:,22)=(MTbef(:,29)+MTbef(:,32))/2;
       MTb(:,23)=(MTbef(:,28)+MTbef(:,33))/2;
       MTb(:,24)=MTbef(:,34);
       
       clear MTbef;
       
       Tb=MTb(1,:);
       Tmr=MTb(2,:);
       Tau=MTb(3,:);
       
return

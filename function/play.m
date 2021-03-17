% Freeware Tic Tac Toe for Matlab 5
% Version 1.0
%
% To begin, simply type: play
% To restart, choose GAME + NEW from the menu (or retype: play)
%
% Tested for Matlab 5 on SGI and PC Windows 95 platforms
% (Version 4 is not supported: e.g. because of switch/case statements)
%
% Program randomly chooses between [geometrically] equivalent strategies
% Sorry, no way to win against the program/computer
%
% By Ivo Penzar, June 1997
% http://www.geocities.com/Heartland/1302/penzar.html
% mailto:penzar@mayo.edu


function play(i,j) 

 t='Ivo''s Tic Tac Toe'; t1='Tic Tac Toe';

 tmp1=get(0,'ShowHiddenHandles');
 set(0,'ShowHiddenHandles','on');
 hf=findobj('Tag',t);
 set(0,'ShowHiddenHandles',tmp1);

 if isempty(hf)
  ver=version;
  k=0; P=[1 0;0 1];
  hf=figure('Tag',t,'Name',t,'NumberTitle','off',...
            'HandleVisibility','callback','IntegerHandle','off');
  hi=uimenu(hf,'Label','&Info');
  hc=uimenu(hi,'Label','&Computer');
  uimenu(hc,'Label',computer);
  if isunix
   [tmp2,lh]=unix('hostname');
   uimenu(hc,'Label',lh);
  end
  hm=uimenu(hi,'Label','&Matlab');
  uimenu(hm,'Label',ver);
  uimenu(hm,'Label',license);
  ha=uimenu(hi,'Label','&Game');
  uimenu(ha,'Label',t1);
  uimenu(ha,'Label','Version 1.0');
  uimenu(ha,'Label','By Ivo Penzar, 1997');
  hg=uimenu(hf,'Label','&Game');
  uimenu(hg,'Label','&New','CallBack','play;');
  zx=0.1; dx=(1-2*zx)/3;
  zy=0.1; dy=(1-2*zy)/3;
  dx2=dx/2; dy2=dy/2;
  x=zx+[0 dx 2*dx]; y=zy+[2*dy dy 0];
  c='x=get(gco,''Userdata'');play(x(1),x(2));'; 
  for i=1:3 
   for j=1:3 
    ho(i,j)=uicontrol(hf,'Units','normalized',...
     'Position',eval(sprintf...
     ('[%f %f %f %f]',x(j),y(i),dx,dy)),... 
     'Userdata',[i,j],'Callback',c); 
   end 
  end 
  u.ho=ho;

 else
  u=get(hf,'Userdata'); 
  k=u.k; P=u.P;
  ho=u.ho; 
 end

 if nargin<2
  k=0; P=[1 0;0 1];
  set(ho,'String',''); 
  if quest(0,t1,hf) 
   k=1;
   P=make(P,random({'','i','j','ij'}));
   [i,j]=perm(P',1,1);
   set(ho(i,j),'String','O'); 
  end

 else 
  ho1=ho(i,j);
  if (k<0)|(~isempty(get(ho1,'String')))
   return
  end
  set(ho1,'String','X');
  [k,P,i,j]=me(k,P,i,j); 
  set(ho(i,j),'String','O'); 
  if (k<0)&quest(-k,t1,hf) 
   return
  end
 end 

 u.k=k; u.P=P; 
 set(hf,'Userdata',u); 
 drawnow; 
 
return 


function [k,P,i,j]=me(k,P,i,j)

 [i,j]=perm(P,i,j); 
 [k,P,i,j]=me1(k,P,i,j); 
 [i,j]=perm(P',i,j);
 
return


function [k,P,i,j]=me1(k,P,i,j)

 switch k

 case 0
  if i==3
   [P,i,j]=make(P,'i',i,j);
  end
  if j==3
   [P,i,j]=make(P,'j',i,j);
  end
  if i>j
   [P,i,j]=make(P,'s',i,j);
  end
  if (i==2)&(j==2)
   k=5;
   P=make(P,random({'','i','j','ij'}));
   i=1; j=1;
  else
   k=2+j;
   i=2; j=2;
  end

 case 1
  if i>j
   [P,i,j]=make(P,'s',i,j);
  end
  if (i==2)&(j==2)
   k=11;
   P=make(P,random({'','s'}));
   i=1; j=3;  
  else
   k=12;
   switch i      
   case 1
    P=make(P,'js');
    j=1; 
   case 2
    i=1;  
   case 3
    P=make(P,random({'','s'}));
    i=1;  
   end
  end

 case 3
  if i>j
   [P,i,j]=make(P,'s',i,j);
  end
  switch i
  case 3
   k=31;
   P=make(P,random({'','s','ij','sij'}));
   i=1; j=2;
  case 2
   k=34;
   i=1; j=3;
  case 1
   k=30+j;
   j=5-j;
  end

 case 4
  if j==3
   [P,i,j]=make(P,'j',i,j);
  end
  switch i
  case 1
   [k,P,i,j]=me1(3,P,1,2);
  case 2
   k=41;
   i=1;
  case 3
   if j==1
    [P,i,j]=make(P,'i',i,j);
    [k,P,i,j]=me1(3,P,3,2);
   else
    k=42;
    P=make(P,random({'','i','j','ij'}));
    i=1; j=1;
   end
  end
 
 case 5  
  if i>j
   [P,i,j]=make(P,'s',i,j);
  end
  switch i
  case 3
   k=51;
   P=make(P,random({'','s'}));
   i=1; j=3;
  case 2
   [P,i,j]=make(P,'s',i,j);
   k=53;
   i=1;
  case 1
   if j==2
    k=52; 
    i=3;  
   else
    k=51;
    [P,i,j]=make(P,'s',3,1);
   end
  end
 
 case 12
  if i==1
   k=121;
   i=3; j=1;
  else
   k=-2;
   i=1; j=2;
  end

 case 11
  if i==1
   k=111;
   i=3;
  else
   k=-2;
   i=1; j=2;
  end

 case 31
  if j==2
   k=311;
   j=1;
  else
   k=-2;
   i=3; j=2;
  end

 case 32
  if (i==3)&(j==1)
   k=321;
   i=2;  
  else
   k=-2;
   i=3; j=1;
  end
 
 case 33
  if j==2
   k=331;
   P=make(P,random({'','j'}));
   i=2; j=3;
  else
   k=-2;
   i=3; j=2;
  end

 case 34
  if (i==3)&(j==1)
   k=341;
   i=2;
  else
   k=-2;
   i=3; j=1;
  end

 case 41
  if i>j
   [P,i,j]=make(P,'s',i,j);
  end
  if i==3
   k=411;
   P=make(P,random({'','s'}));
   i=1;
  else
   k=-2;
   i=3;
  end

 case 42
  if (i==3)&(j==3)
   k=421;
   j=1;
  else
   k=-2;
   i=3; j=3;
  end

 case 51
  if i==1
   k=511;
   i=3;
  else
   k=-2;
   i=1; j=2;
  end

 case 52
  switch i
  case 1
   k=522;
   i=3; j=1;  
  case 2
   k=520+j;
   j=4-j;  
  case 3
   if j==1
    k=511;     
    [P,i,j]=make(P,'j',i,j);
   else
    k=531;
    [P,i,j]=make(P,'i',i,j);
   end
   i=1; j=1;
  end

 case 53
  if i==1
   k=531;
   i=3; j=1;
  else
   k=-2;
   i=1; j=3;
  end

 case 121
  k=-2;
  i=2;  
  if j==1
   j=2;
  else
   j=1;
  end
 
 case 111
  if j==1
   [P,i,j]=make(P,'j',i,j);
  end
  k=1109+i;
  i=2; j=1;
   
 case 311
  if i==1
   k=-1;
  else
   k=-2;
   j=3;
  end
  i=3-i;
  
 case 321
  if i==2
   k=-1;
  else
   k=-2;
   j=3;
  end
  i=5-i;
  
 case 331
  if i==2
   k=-1;
  else
   k=-2;
   j=1;
  end
  i=5-i;
  
 case 341
  k=-1;  
  if i==3
   j=5-j;
  else
   i=3;
  end
  
 case 411
  if j==1
   k=-1;
   j=2;  
  else
   k=-2;
   i=3; j=1;
  end
   
 case 421
  k=-2;  
  if j==1
   i=1; j=3;
  else
   i=2; j=1;
  end
 
 case 511
  k=-1;  
  if i==2
   j=4-j;
  else
   P=make(P,random({'','j'}));
   i=2;
  end

 case 521
  if i>j
   [P,i,j]=make(P,'s',i,j);
  elseif i==3
   P=make(P,random({'','s'}));
  end
  k=-1;
  i=3; j=1;
  
 case 522
  k=-2;  
  if i==2
   i=3; j=3;
  else
   i=2; j=1;
  end
 
 case 523
  if j==1
   k=-1;
   i=1;
  else
   k=-2;
   i=3;
  end
  j=4-j;
  
 case 531
  if j==1
   k=-1;
  else
   k=-2;
   i=2;
  end
  j=4-j;

 case {1111,1112}
  if j==1
   i=1114-k;  
   k=-1;
  else
   k=-2;
   i=3;
  end
  j=4-j;
  
 end

return


function [P,i,j]=make(P,s,i,j)

 l=length(s);

 n=l;
 while n>0
  switch s(n)
  case 'i'
   P(1,:)=-P(1,:);
  case 'j'
   P(2,:)=-P(2,:);
  case 's'
   tmp1=P(1,:);
   P(1,:)=P(2,:); P(2,:)=tmp1;
  end
  n=n-1;
 end

 if nargin>=4
  n=l;
  while n>0
   switch s(n)
   case 'i'
    i=4-i;
   case 'j'
    j=4-j;
   case 's'
    tmp2=i;
    i=j; j=tmp2;
   end
   n=n-1;
  end
 end

return


function [i,j]=perm(P,i,j)

 tmp=P*[i-2;j-2];
 i=tmp(1)+2; j=tmp(2)+2;

return


function e=random(x)

 e=x{ceil(length(x)*rand)};

return


function rtn=quest(k,t,h)

 if k
  q={'Let''s play again...','Sorry, you''ve lost...'};
  o1='New'; o2='Quit'; o3='Cancel';
  switch questdlg(q{k},t,o1,o2,o3,o1);
  case {o1}
   play;
   rtn=1;
  case {o2}
   delete(h);
   rtn=2;
  case {o3}
   rtn=0;
  end
  
 else
  o1='Me'; o2='You';
  switch questdlg('Starting?',t,o1,o2,o2);
  case {o1}
   rtn=1;
  case {o2}
   rtn=0;
  end
 end   
 
return


% dirsort gives back the list of file in 
% a directory, sorted depending on user choice.
%
% Es:
%    [outlist]=dirsort(where,how)
% Input: 
%    where: directory path (wildcard for file names are accepted)
%    how:   how to sort files ('name','size','date')
% Output:
%    outlist: a structure conteining fields: name, date, bytes and isdir
%
% Nico, 10/2000

function [outlist]=dirsort(where,how)

  listofile=dir(where);  
    
  switch how
    case 'name'
      [outlist]=sortname(listofile); 
    case 'size'
      [outlist]=sortsize(listofile); 
    case 'date'
      [outlist]=sortdate(listofile); 
    otherwise
      fprintf(1,'What does mean "%s" ? Sorry, abort sorting.\n',how);
  end
  
  
return

%%%% FUNCTIONS %%%%

function [outlist]=sortname(inlist)  

  hmany=length(inlist);
  acell=cell(hmany,1);

  for nf=1:hmany     
     acell(nf)={inlist(nf).name};
  end   
  
  [junk,indx]=sort(acell);
  outlist=inlist(indx);

return

%
function [outlist]=sortsize(inlist)  

  hmany=length(inlist);
  bytes=zeros(hmany,1);

  for nf=1:hmany     
     bytes(nf)=inlist(nf).bytes;
  end   
  
  [junk,indx]=sort(bytes);
  outlist=inlist(indx);

return

%
function [outlist]=sortdate(inlist)  

  hmany=length(inlist);
  stime=zeros(hmany,1);

  for nf=1:hmany     
     stime(nf)=datenum(inlist(nf).date);
  end   
  
  [junk,indx]=sort(stime);
  outlist=inlist(indx);

return

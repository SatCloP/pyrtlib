% THIS IS A ROUTINE FOR COUNTING NUMBER OF ROWS IN A RAOB
% YOU HAVE TO GIVE IN INPUT 
% fid    file id
% nskip  number of number or string to skip in starting up
% nc     number of raob's colons 
% flg    checking flag that means end of raob (string or number)
% nbreak colon position in a row of checking flag 
function [nl]=counts(fid,nskip,nc,flg,nbreak)

sorn=class(flg);
if length(sorn)<6
   sorn(1,5:6)='  ';
end

if sorn=='double'
   str=0;
   item=0;
elseif sorn=='char  '
   str=1;
   item(1,length(flg))=' ';
end

nl=0;

skip=fscanf(fid,'%s',[1 nskip]);
while item~=flg
   nl=nl+1;
   try
      MTX(nl,:)=fscanf(fid,'%f',[1 nc]);
      item=MTX(nl,nbreak);
      if str
         item=num2str(item);
         if length(item)>length(flg)
            item=item(1,1:length(flg));
         elseif length(item)<length(flg)
            item(1,length(item)+1:length(flg))=' ';
         end
      end   
   catch
      item=flg;
   end
end

nl=nl-1;

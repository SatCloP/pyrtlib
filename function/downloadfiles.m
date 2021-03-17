% This function was written for WViop2000 and allows 
% dowloading data from remote computer to local.
% I will improve it to be more general
%

function downloadfiles

  %kind=['2ch'; 'met'; 'psr'; '5mm'; 'enc'];
  %kind=['2ch'; 'met'; 'psr'; 'enc'];
  kind=['enc';];

    !ftp -vs:ftp_filenames 140.172.33.1
  
  for ik=1:length(kind(:,1))
     
     list2down=matchfiles(kind(ik,:));
     if isempty(list2down); 
        fprintf('No new %s file to download...\n',kind(ik,:)); 
     else   
        fprintf('%4d %s file to download...\n',length(list2down(:,1)),kind(ik,:)); 
        writeftp(list2down,kind(ik,:));
        cd([kind(ik,:) '/data'])     
         %! ftp -vs:ftp_file 140.172.33.1
        cd ../../
     end;
pause
  end
  
return


%%%%%%%%%% Functions %%%%%%%%%%%%%

function list2down=matchfiles(what)

    list2down=[];
    fid=fopen([what '_fnames.dat'],'r');
    listoftaken=dir([what '\data']);
    taken=strvcat(listoftaken.name);
    while ~feof(fid)
       fname=fscanf(fid,'%s',1);
       %ismember(fname,taken,'rows');
        %if feof(fid); return; end; 
        if ~ismember(fname,taken,'rows')
           list2down=strvcat(list2down,fname);
        end
    end
    fclose(fid);
return

%
function writeftp(list,what)

  str=['/hd1/ok00/' what '/data'];
  
  fid=fopen('ftp_file','w');
  fprintf(fid,'csr\ncsr\nprompt\nbin\ncd %s\n',str);
  for in=1:length(list(:,1))
    fprintf(fid,'get %s\n',list(in,:));
  end  
  fprintf(fid,'bye');
  fclose(fid);
    
return

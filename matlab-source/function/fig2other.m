% fig2other converts Matlab figure filename.fig
% in a different format with the same name.
% It isn't a smart idea but save time
% if the number of file is big.
%
% N.B.: Do not put .fig at the end of wildcard
%
% Usage:
%    fig2other(dir,extention,wildcard)
% Es:
%    fig2other('C:\Nico\Matools\','jpg')
%    fig2other(pwd,'jpg','sa*')

function fig2other(where,how,wild)

if where(end)~='/'; where=[where '\']; end;
if nargin < 3; wild=''; end;

listofig=dirsort([where wild '*.fig'],'name');
numfig=length(listofig);

if ~numfig; fprintf(1,'No .fig file in %s \n',where); return; end;

for nf=1:numfig
   
   jnk=findstr(listofig(nf).name,'.');
   figname=listofig(nf).name(1:jnk(end)-1); % there might be more than one "."
   fprintf(1,'Converting figure %s .....',figname);
   
   open([where figname '.fig']);
   format4paper(gcf,11);
   pause(1);
   saveas(gcf,[where figname '.' how],how);
   pause(1);
   close(gcf);
   fprintf(1,' Ok.\n');
      
end

fprintf(1,'Converted successfully %d figure(s)\n',numfig);

return
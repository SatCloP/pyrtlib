% FUNCTION THAT ADDS A LINE TO THE LEGEND OF THE CURRENT FIGURE.
% X ES.: addlegend(newstrg)
%
% IF IT GIVES ERRORS TRY TO CHANGE "if 1" IN "if 0"
% I COULD HAVE TESTED JUST FEW TIMES
%
% Nico,2000
%
% Hystory:     08/04/2000
% Last update: 2004/12/13 - changed to allow adding to 2 or more strings


function [lstrg]=addlegend(newstrg,legendpos)

lstrg=[]; % avoiding "Warning: One or more output arguments not assigned.."

if nargin < 2
   legendpos=1;
end 

if 1
  hc=get(gcf,'Children');
  lstruct=get(hc(1),'Userdata');
  if isempty(lstruct);
     if ~isempty(newstrg)
        legend(newstrg); 
     end  
     return; 
  end;
  lstrg=lstruct.lstrings;
  lstrg(end+1)={newstrg};
  %lstrg=get(lstruct.LabelHandles(1),'String');
  %try
  %  lstrg=strvcat(lstrg,newstrg);
  %catch
  %  lstrg(end+1)={newstrg};
  %end  
  legend off
  legend(lstrg,legendpos);
else
  hl=findobj('Tag','legend');
  lstruct=get(hl,'Userdata');
  if isempty(lstruct); return; end;
  lstrg=get(lstruct.LabelHandles(1),'String');
  lstrg=strvcat(lstrg,newstrg);
  legend off
  legend(lstrg,legendpos);
end
  
return
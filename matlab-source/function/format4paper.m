% THIS FUNCTION TRANSFORM A MATLAB DEFAULT FIGURE
% IN A DIFFERENT FORMAT, MORE APPROPRIATE FOR
% PAPERS AND VIEWGRAPHS
% IT NEEDS IN INPUT A FIGURE HANDLE
%
% format4paper(hfig,fsize,msize,lwidth) 
%
% 2005/12/12 - Modifyed to use strcmp instead of if str=='whatever'
% 2006/01/09 - Modifyed to accept multiple figure handles
% 2006/02/23 - Modifyed to accept colormap with colorbars
% 2006/04/27 - Modifyed to accept markersize
% 2009/01/15 - Modifyed to debug
% 2011/01/21 - Modifyed to accept line width
% 2011/01/28 - Modifyed to work with double y axes
% 2013/05/28 - Modifyed to extend markersize to all types of markers

function format4paper(hfig,fsize,msize,lwidth)

switch nargin 
    case 1
      fsize = 14;
      msize = 10;
      lwidth = 1.5;
    case 2
      msize=10;
      lwidth = 1.5;
    case 3
      lwidth = 1.5;
end

for nf = 1:length(hfig)
    
figure(hfig(nf))

plotbox=[];
legends=[];
clrbars=[];
hc=get(gcf,'Children');
for ic = 1:length(hc)
   %if isempty(get(hc(ic),'Tag')) % 2009/01/15
   if isempty(get(hc(ic),'Tag')) & strcmp(get(hc(ic),'Visible'),'on')
      plotbox=[plotbox ic];    
   elseif strcmp(get(hc(ic),'Tag'),'legend');
      legends=[legends ic];
   elseif strcmp(get(hc(ic),'Tag'),'Colorbar');
      clrbars=[clrbars ic];
   end
end

% plots
for iax = 1:length(plotbox)

%   axes(hc(plotbox(iax)));
%   hc1=get(gca,'Children');
   ca1=hc(plotbox(iax));
   hc1=get(ca1,'Children');
   
   for ic1 = 1:length(hc1)
       if strcmp(get(hc1(ic1),'Type'),'line')
          set(hc1(ic1),'Linewidth',lwidth);
       elseif strcmp(get(hc1(ic1),'Type'),'text')
          %pos=get(hc1(ic1),'position'); pos(1)=pos(1)-0.5;
          %set(hc1(ic1),'position',pos);
          set(hc1(ic1),'Fontsize',fsize,'Fontweight','bold');
       end
   end
    
   set(ca1,'Fontsize',fsize,'Fontweight','bold');

   htl=get(ca1,'Title');
   set(htl,'Fontsize',fsize,'Fontweight','bold')
   hxl=get(ca1,'Xlabel');
   set(hxl,'Fontsize',fsize,'Fontweight','bold')
   hyl=get(ca1,'Ylabel');
   set(hyl,'Fontsize',fsize,'Fontweight','bold')

end

% legends
%hc=get(gcf,'Children');
for il=1:length(legends)
   axes(hc(plotbox(il)));
   lstruct=get(hc(legends(il)),'Userdata');
   if ~isempty(lstruct);
      legend off;
      legend(lstruct.lstrings,lstruct.legendpos)
   end
end

% colorbars
if 0 % it does not work on new matlab versions
for il=1:length(clrbars)
   axes(hc(clrbars(il)));
   set(gca,'Fontsize',fsize,'Fontweight','bold');
   htl=get(gca,'Title');
   set(htl,'Fontsize',fsize,'Fontweight','bold');
   hxl=get(gca,'Xlabel');
   set(hxl,'Fontsize',fsize,'Fontweight','bold');
   hyl=get(gca,'Ylabel');
   set(hyl,'Fontsize',fsize,'Fontweight','bold');
end
end

% markersize
set(findobj('Marker','.'),'Markersize',msize);
set(findobj('Marker','+'),'Markersize',msize);
set(findobj('Marker','x'),'Markersize',msize);
set(findobj('Marker','o'),'Markersize',msize);
set(findobj('Marker','*'),'Markersize',msize);
set(findobj('Marker','^'),'Markersize',msize);
set(findobj('Marker','v'),'Markersize',msize);
set(findobj('Marker','s'),'Markersize',msize);
set(findobj('Marker','d'),'Markersize',msize);

end % end of loop on figure handles

return
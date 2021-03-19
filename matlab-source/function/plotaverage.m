% FUNCTION THAT AVERAGES THE CURVES IN A PLOT,
% GIVING BACK THE SAME PLOT, BUT WITH AVERAGED POINTS
%
% X Es: plotaverage(np,add,usrmark)
% 
% IT NEEDS THE NUMBER OF POINTS TO AVERAGE (np), A FLAG FOR
% ADDING NEXT PLOT OR NOT ('Y' OR 'N') AND THE MARKER STYLE
% (OR JUST '' TO USE THE ORIGINAL ONES).
%
% Nico, 2000


function plotaverage(np,add,usrmark)

  hc=get(gca,'Children');
  lgndstr=addlegend([]);
  
  for ic=1:length(hc)
     
     X=get(hc(ic),'Xdata');
     Y=get(hc(ic),'Ydata');
     clr=get(hc(ic),'Color');
     mrk=get(hc(ic),'Marker');
     mrksz=get(hc(ic),'Markersize');
     
     XA=boxmean(X',np);
     YA=boxmean(Y',np);
     
     FEATs(ic)=struct('Xdata',X,'Ydata',Y,'XA',XA,'YA',YA,'Color',clr,'Marker',mrk,'Markersize',mrksz);

  end
  
  if ~strcmp(usrmark,''); [FEATs.Marker]=deal(usrmark); end;
  if strcmp(add,'N')
     figure
  end
  hold on
    
  for ic=1:length(hc)
     
    icr=length(hc)-ic+1; % icreverse: it's necessary to match legend labels
    plot(FEATs(icr).XA,FEATs(icr).YA,'Color',FEATs(icr).Color,'Marker',FEATs(icr).Marker,'Markersize',FEATs(icr).Markersize)

  end
  box on
  zoom on
  
  if ~isempty(lgndstr)
     addlgndstr=strcat(lgndstr,[' (' num2str(np) ' points average)']);
     if strcmp(add,'Y'); addlgndstr=strvcat(lgndstr,addlgndstr); end;    
     legend(addlgndstr);
  end
  
return
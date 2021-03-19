function [ax1,ax2,hl1,hl2] = plotxx(x1,y1,x2,y2);
%PLOTXX  Plots two properties as a function of Y 
%
%Similar to PLOTYY, but ...
%the independent variable is on the y-axis, 
%and both dependent variables are on the x-axis.
%
%Syntax: [ax1,ax2,hl1,hl2] = plotxx(x1,y1,x2,y2);
%
%Inputs:  X1,Y1 are the data for the first line (black)
%         X2,Y2 are the data for the second line (red)
%
%The optional output handle graphics objects AX1,AX2,HL1,HL2
%allow the user to easily change the properties of the plot.
%
%Example: Plot temperature T and salinity S 
%         as a function of depth D in the ocean
%
% [axT,axS,hlT,hlS] = plotxx(T,D,S,D);


%The code is inspired from page 10-26 (Multiaxis axes)
%of the manual USING MATLAB GRAPHICS, version 5.
%
%Tested with Matlab 5.2.1.1420 on PCWIN

%Author: Denis Gilbert, Ph.D., physical oceanography
%Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%November 1997; Last revision: 18-JAN-1999

hl1=line(x1,y1,'Color','k');
ax1=gca;
set(ax1,'Position',[0.15 0.15 0.75 0.70])
set(ax1,'XColor','k','YColor','k');

ax2=axes('Position',get(ax1,'Position'),...
   'XAxisLocation','top',...
   'YAxisLocation','right',...
   'Color','none',...
   'XColor','r','YColor','k');

hl2=line(x2,y2,'Color','r','Parent',ax2);

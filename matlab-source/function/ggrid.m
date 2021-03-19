function x=ggrid(delx,dely)
%	ggrid.m draws in the major grid lines only.
%
%	Syntax:
%	x=ggrid(delx,dely)
%
%   will draw grid lines at certain values.  Just type
%   "ggrid" w/o arguments to
%   draw lines at the decades, which is useful on a
%   log-log plot.
%
%   Enter arguments
%   delx to control x major grid spacing,
%   dely to control y major grid spacing.
%
%   Output is the plot extremes: a = [xmin,xmax,ymin,ymax]
%
%	Example:
%	ggrid - to put vertical & horzontal lines every decade.
%	ggrid(500,.02) put vertical lines every 500 units and
%                        horizontal ever .02 units.
%
%	Try typing in "grid" then, after the grid is displayed,
%	type "ggrid".
%
%
%	Uses the get(gca,'Xcolor') to determine what
%	color lines to use.
%
%   Bill Dunn 9/96 org 9735
%   wndunn@sandia.gov
%
%  SEE ALSO GRID, LINE
%
lcolor=get(gca,'Xcolor');
%
% if no arguments in, use decade spacing.
	a=axis;
if nargin == 0
	logrange(1) = floor(log10(a(1)));
	logrange(3) = floor(log10(a(3)));
	logrange(2) = ceil(log10(a(2)));
	logrange(4) = ceil(log10(a(4)));
	% draw vertical lines
	for i=logrange(1):logrange(2)
		line('Color',lcolor,'XData',[10^i;10^i],'YData',[a(3);a(4)])
	end
	% draw horizontal lines
	for i=logrange(3):logrange(4)
		line('Color',lcolor,'XData',[a(1);a(2)],'YData',[10^i;10^i])
	end
	line('Color',lcolor,'XData',[a(2);a(2)],'YData',[a(3);a(4)])	
	line('Color',lcolor,'XData',[a(1);a(2)],'YData',[a(4);a(4)])
end
% if 2 arguments (xspacing & yspacing) are entered:
if nargin == 2
	% draw vertical lines
	for i=a(1):delx:a(2)
		line('Color',lcolor,'XData',[i;i],'YData',[a(3);a(4)])
	end
	% draw horizontal lines
	for i=a(3):dely:a(4)
		line('Color',lcolor,'XData',[a(1);a(2)],'YData',[i;i])
	end
end

if nargin == 1
	disp(['Need to enter 0 or two arguments'])
end
return


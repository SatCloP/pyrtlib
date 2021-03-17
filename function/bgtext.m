function handle = bgtext(x, y, z, string, color, varargin)

% bgtext -- Text with background color.
%  bgtext(x, y, 'string', color, ...) draws the text with
%   the given patch color in the background.  The "ResizeFcn"
%   of the figure is set to resize the patches automatically.
%   The default color is light-gray.  Except for the added
%   "color" argument, the syntax is the same as for "text".
%  bgtext(x, y, z, 'string', color, ...) for 3-d text position.
%  bgtext(theFigure) resizes the background patches to fit the
%   extents of the "bgtext" objects in theFigure (default = gcf).
%  bgtext('demo', N) demonstrates itself for N random locations
%   (default = 2), in a random color.
 
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
 
% Version of 26-Aug-1999 14:26:10.
% Updated    26-Aug-1999 16:46:22.

DEFAULT_COLOR = [9 9 9]/10;

if nargin < 1, x = gcf; end

% Demonstration.

if isequal(x, 'demo')
	n = 2;
	if nargin > 1
		if isstr(y), y = eval(y); end
		n = y;
	end
	x = rand(n, 1);
	y = rand(size(x));
	string = 'Hello World';
	color = 0.5 * (1 + rand(1, 3));
	color = color / max(color);
	set(gca, 'XLim', [-0.5 1.5], 'YLim', [-0.5 1.5])
	bgtext(x, y, string, color)
	set(gcf, 'Name', [mfilename ' demo'])
	zoomsafe
	figure(gcf)
	return
end

% Resize background patch sizes.

if nargin < 2 & ishandle(x) & isequal(get(x, 'Type'), 'figure')
	h = x;
	f = findobj(h, 'Type', 'text', 'Tag', mfilename);
	for i = 1:length(f)
		t = f(i);
		p = get(t, 'UserData');
		if ishandle(p)
			switch get(p, 'Type')
			case 'patch'
				extent = get(t, 'Extent');
				x0 = extent(1); y0 = extent(2);
				width = extent(3); height = extent(4);
				xp = x0 + [0 width width 0];
				yp = y0 + [0 0 height height];
				set(p, 'XData', xp, 'YData', yp)
			otherwise
			end
		end
	end
	return
end

% Normal processing.

is3d = ~isstr(z);

if ~is3d
	if nargin > 4
		varargin = [{color} varargin];
	else
		varargin = {};
	end
	if nargin > 3
		color = string;
	else
		color = DEFAULT_COLOR;
	end
	string = z;
else
	if nargin < 5
		color = DEFAULT_COLOR;
	end
end

wasHold = ishold;
hold on

% Measure temporary text (invisible).

if ~is3d
	t = text(x, y, string);
else
	t = text(x, y, z, string);
end

if ~isempty(varargin), set(t, varargin{:}), end
set(t, 'Visible', 'off')

extent = zeros(length(t), 4);
for i = 1:length(t)
	extent(i, :) = get(t(i), 'Extent');
end

% Delete temporary text.

delete(t)

% Create patches.

p = zeros(length(t), 1);
for i = 1:length(t)
	x0 = extent(i, 1); y0 = extent(i, 2);
	width = extent(i, 3); height = extent(i, 4);
	xp = x0 + [0 width width 0];
	yp = y0 + [0 0 height height];
	p(i) = patch(xp, yp, color);
end
set(p, 'EdgeColor', 'none', 'Tag', mfilename)

% Superimpose new text.

if ~is3d
	t = text(x, y, string);
else
	t = text(x, y, z, string);
end

if ~wasHold, hold off, end

if ~isempty(varargin), set(t, varargin{:}), end

for i = 1:length(t)
	set(t(i), 'UserData', p(i), 'Tag', mfilename)
end

% Set the "ResizeFcn" callback.

set(gcf, 'ResizeFcn', [mfilename])

% Return.

if nargout > 0, handle = t; end

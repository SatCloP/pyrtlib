%function stackfigs
% <cpp> stack figure windows


%Charles Plum                    Nichols Research Corp.
%<cplum@nichols.com>             70 Westview Street
%Tel: (781) 862-9400             Kilnbrook IV
%Fax: (781) 862-9485             Lexington, MA 02173


maxpos  = get (0,'screensize'); % determine terminal size in pixels
hands   = get (0,'Children');   % locate all open figure handles
hands   = sort(hands);          % sort figure handles
numfigs = size(hands,1);        % number of open figures
maxfigs = fix(maxpos(4)/20);


if (numfigs>maxfigs)            % figure limit check
        disp([' More than ' num2str(maxfigs) ' requested '])
        return
end


% tile figures by postiion 
% Location (1,1) is at bottom left corner
pnum=0;
for iy = 1:numfigs
  figure(hands(iy))
  p = get(hands(iy),'Position'); % get figure position
  
  ypos = maxpos(4) - (iy-1)*20 -p(4) -45 ; % figure location (row)
  xpos = fix((iy-1)*5 + 15);     % figure location (column)
  
  set(hands(iy),'Position',[ xpos ypos p(3) p(4) ]); % move figure
  
end
return
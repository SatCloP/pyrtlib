% textp places the text proportionally to the x and y limits
%
% Usage:
%        [xt,yt,ht]=textp(ha,xp,yp,strng);
%
% Inputs:
%        ha: axis handle 
%        xp: fraction of the x axis where to start the text (es: 2/3)
%        xp: fraction of the y axis where to start the text (es: 1/3)
%      strg: string to use with text function
% Outputs:
%        xt: x value where text started
%        yt: y value where text started
%        ht: text handle (for changing the color, font, etc..)
% Es:
%        [xt,yt,ht]=textp(gca,2/3,1/3,'CIAO');
%
% Nico, Oct 2001

function [xt,yt,ht]=textp(ha,xp,yp,strng)

xl = get(ha,'xlim');
yl = get(ha,'ylim');

xt = xl(1) + xp * ( xl(2)-xl(1) );
yt = yl(1) + yp * ( yl(2)-yl(1) );

ht = text(xt,yt,strng);

return
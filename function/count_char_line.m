% count_char_line(n) prints on screen two strings that may be useful for
% counting characters in a line.
%
% Ex:
%    count_char_line(3) prints out:
%000000000111111111122222222223
%123456789012345678901234567890
%
% Nico, Mar 2010
%
% History:
% 2012/12/14 - Modified to fix a bug for n>9
% 2014/08/04 - Completely rewritten to print 3 lines for n>9 

function count_char_line(n)

if nargin < 1
   n = 9;
end

% for in = 1:n
%     while in > 9; in = in-10; end;
%     fprintf(1,'---------%1d',in);
% end
% fprintf(1,'\n');
% 
% for in = 1:n
%     fprintf(1,'1234567890');
% end
% fprintf(1,'\n');
   
MTX = 1:n*10;
nc = 2; if n > 9; nc = 3; end;

for in = 1:length(MTX)
    STR(in,1:nc) = twodigstr(MTX(in),nc);
end

STR'

return
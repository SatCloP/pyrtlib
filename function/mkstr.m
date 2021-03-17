% mkstr
% This function initializates a matrix of empty strings given in input the 
% number of rows and columns
% Es: str = mkstr(5,10);
%
% Nico, Apr 2007

function outstr = mkstr(nr,nc)

for ir = 1:nr
    for ic = 1:nc
        outstr(ir,ic) = ' ';
    end
end

end
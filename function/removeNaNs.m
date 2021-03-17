function newAP = removeNaNs(ALLPLANES)
% REMOVENANS - removes NaN rows and columns from data which has been interpolated.
%
% a = [ NaN NaN NaN     removeNaNs(a) will return:  newAP = [ 1 2
%       1   2   NaN                                           3 4
%       3   4   NaN                                           5 6
%       5   6   NaN                                           7 8 ]
%       7   8   NaN ]

% Author: Andrew Hastings
% Created: Feb. 19, 1999
% Email: andrewh@qualcomm.com
% Tested with Matlab v5

[i,j] = find(isfinite(ALLPLANES));		% Get the indicies of the actual data.
i = unique(i);									% Need to run unique on them as the values are repeated.
j = unique(j);								
endMatrix = []; idx = 1;					% Initialize endMatrix and index.

if ~isempty(find(isfinite(ALLPLANES)))
  pages = j(diff(j)~=1);					% The number of z planes-1
  for numPg=1:length(pages)+1			
   
    if isempty(pages)						% If there's only 1 page
      cols = length(j);						% then the number of columns per page = length of j
    else 
      cols = pages(1);						% Otherwise it'll be equal to the first number in pages.
    end
  
    pagesAt=[idx:cols*numPg];				% Get the next page (set of z plane data) of j indicies.
    idx = cols*numPg+1;					
    endMatrix=cat(3,endMatrix,ALLPLANES(i,j(pagesAt)));
  end
end
newAP=endMatrix;
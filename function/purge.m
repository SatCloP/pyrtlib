% GIVING IN INPUT A MATRIX, A FLAG AND A COLOUMN 
% NUMBER THIS FUNCTION PURGES THE MATRIX LINES IN WHICH 
% MATRIX(NLINE,NCOLOUMN) IS EQUAL TO THE FLAG VALUE 
% Es:
%    [MTXP]=purge(MTX,nc,flag);

function [MTXP]=purge(MTX,nc,flag);

MTXP=MTX;

ip=find(MTX(:,nc)==flag);

if ~isempty(ip)
  MTXP(ip,:)=[];
end

return
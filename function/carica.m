% Function that allows to load matrix from a file,
% forgetting the unconfortables features of fscanf.
%
% Es:
%    MX=carica(fid,nc,nr,format);
% Inputs:
%    fid: file id
%    nc : number of colons (as needed from the fid file)
%    nr : number of rows   (as needed from the fid file)
%    format: usual Matlab string for format ('%f','%s'...)

function MX=carica(fid,nc,nr,format);

if nr==0
   MX=(fscanf(fid,format,[nc,inf]))';   
else   
   MX=(fscanf(fid,format,[nc,nr]))';
end

return





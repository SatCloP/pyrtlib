% This function takes in input two structures with the same fields
% and gives back one structure with the same fields filled with
% fields values of the first and the second structure.
%
% N.B.: Right now it works only if the fields are row vectors (as farran matlab files)
% N.B.: It should work with both row and col concatenation
%
% Es:
%    [tot]=meltstructs(struct1,struct2); 
%
%    [tot]=meltstructs(struct1,struct2,'col'); 
%
% Nico, Jun 2001
% 2010/05/07 - Changed to work with 'row' and 'col' concatenation

function [tot]=meltstructs(struct1,struct2,how)

if nargin < 3
    how = 'row'
end

tot=struct1;
fnames = fieldnames(struct1);

for nf=1:length(fnames)
   
   fname=char(fnames(nf));
   
   field1=getfield(struct1,fname);
   field2=getfield(struct2,fname);
   tot = setfield(tot,fname,[]);
   switch how
       case 'row'
          totfield = [field1 field2];
       case 'col'
          totfield = [field1; field2];
   end 
   tot = setfield(tot,fname,totfield);
   
end

return
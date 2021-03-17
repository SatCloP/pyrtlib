% FUNCTION TO WRITE A TWO DIMENSIONAL MATRIX IN A FILE
%
% ES: writeinfile(fid,V,form)
%
% fid   file identification
% V     Matrix
% form  format (could be just a string, like ' %f', or a colon vector of strings)

function writeinfile(fid,V,form)

for il=1:length(V(:,1))
   for ic=1:length(V(1,:))
      if length(form(:,1))==1
         fprintf(fid,form,V(il,ic));
      else
         fprintf(fid,form(ic,:),V(il,ic));
      end   
   end
      fprintf(fid,'\n');
end

return
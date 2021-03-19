% MATRIX LINEAR REGRESSION FUNCTION
% UPGRADED WITH EXTERNAL TEST CHOICE
% In input matrix each colon must correspondes to a different 
% variable while each row to a different observation.
% Le matrici in input devono avere le colonne corrispondenti 
% alle varie grandezze e le righe ai vari conteggi
%
% [Vy_r]=linreg(Vx,Vy,vtx);
%
% N.B.: IT'S OBSOLETE NOW. NEW VERSION IS linreg2.m 
%       BUT STILL SOME PROGRAM MIGHT USE THIS ONE.

function [Vy_r]=linreg(Vx,Vy,vtx);

Vxm=mean(Vx);
dimVx=ones(length(Vx(:,1)),1);

Vym=mean(Vy);
dimVy=ones(length(Vy(:,1)),1);

vtxm=mean(vtx);
dimvtx=ones(length(vtx(:,1)),1);

sgmVx=(std(Vx)); 
Vxsc=(Vx-dimVx*Vxm)./(dimVx*sgmVx);

Vysc=Vy-dimVy*Vym;

CVx=cov(Vxsc);
CRCVxy=(Vysc'*Vxsc)/(length(dimVx)-1);

D=CRCVxy*inv(CVx);

if length(vtx)==1 & vtx(1)==0
  Vy_r=dimVy*Vym+Vxsc*D';
else 
  vtxsc=(vtx-dimvtx*Vxm)./(dimvtx*sgmVx);
  Vy_r=dimvtx*Vym+vtxsc*D';
end

return
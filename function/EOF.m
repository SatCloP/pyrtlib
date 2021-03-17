% THIS FUNCTION COMPUTES THE EMPIRICAL ORTHOGONAL FUNCTIONS (EOF).
% INPUTS ARE: MEASURE MATRIX, THRESHOLD FOR EIGENVALUES OF COVARIANCE
% MATRIX AND NOISE LEVEL TO ADD AT MEASURES.
% INPUT MATRIX MUST HAVE ONE VARIABLE FOR EACH COLUMN AND ONE ROW FOR
% EACH MEASURE.
% OUTPUTS ARE: COEFFICIENT MATRIX AND EOF MATRIX
% [V,D] = EIG(X) produces a diagonal matrix D of eigenvalues and a
% full matrix V whose columns are the corresponding eigenvectors so
% that X*V = V*D.
%
% [Cs,eof,perc]=EOF(M,thr,noise);
%
% Last update: 2004/05/11

function [Cs,eof,perc,egvs]=EOF(M,thr,noise);

CV=cov(M);
[Vs,As]=eig(CV);

egvs=diag(As);
[egvs,inxs]=sort(-egvs);
egvs=-egvs;

for i=1:length(egvs)
    perc(i) = sum(egvs(1:i))/sum(egvs);
end

Vs=Vs(:,inxs);

ncp=find(egvs>thr);

eof=Vs(:,ncp);
Mn=M+noise*randn(size(M));
Cs=Mn*eof;



return
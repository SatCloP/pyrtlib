% UPGRATED VERSION (26/04/00) MATRIX LINEAR REGRESSION FUNCTION
% 
% In the Input matrix each colomn must correspondes to a different 
% variable while each row to a different observation.
% (Le matrici in input devono avere le colonne corrispondenti 
% alle varie grandezze e le righe ai vari conteggi)
%
% TWO DIFFERENT USAGES:
% 
% 1) CALCULATING LINEAR REGRESSION COEFFICIENT MATRIX
%    [D,Vxm,sgmVx,Vym]=linreg2(Vx,Vy)
%
% 2) COMPUTING REGRESSION
%    [vy_r]=linreg2(D,Vxm,sgmVx,Vym,vx)
%
% LEGEND:
%
% Vx:    Vectors of Radiometric (and Surface) Measurements
% Vy:    Vectors of Atmospheric Parameters Measurements
% V?m:   Average of V? (Vx or Vy)
% sgmVx: Standard Deviation of Vx 
% D:     Matrix of Retrieval Coefficients
% vx:    Indipendent Measurements (out of DataBase)
% vy_r:  Retrieved Atmospheric Parameters
%
% Nico Cimini, 04/2000


function [varargout]=linreg(varargin);

switch nargin
   
case 2
   
   A=varargin{1};
   B=varargin{2};
   
   Am=mean(A);
   dimA=ones(length(A(:,1)),1);

   Bm=mean(B);
   dimB=ones(length(B(:,1)),1);

   sgmA=(std(A)); 
   Asc=(A-dimA*Am)./(dimA*sgmA);

   Bsc=B-dimB*Bm;

   CV=cov(Asc);
   CRCV=(Bsc'*Asc)/(length(dimA)-1);

   D=CRCV*inv(CV);

   varargout{1}=D;
   varargout{2}=Am;
   varargout{3}=sgmA;
   varargout{4}=Bm;

case 5
   
    D=varargin{1};
    Am=varargin{2};
    sgmA=varargin{3};
    Bm=varargin{4};    
    a=varargin{5};
    
    dima=ones(length(a(:,1)),1);

    asc=(a-dima*Am)./(dima*sgmA);
    a_r=dima*Bm+asc*D';
    
    varargout{1}=a_r;
    
otherwise
   
   fprintf(2,' Wrong number of input variables: %d \n',nargin);
   
end

return
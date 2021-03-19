function [a,b,asig,bsig,chisquared,q]=chi2(x,y,xsig,ysig);

% function [a,b,asig,bsig,chisquared,q]=chi2(x,y,xsig,ysig);
%
% Determines the best fit straight line (a+bx) to the x,y data with
% x and y uncertainties, xsig,ysig.  asig and bsig are the weighted
% uncertainties in a and b.  chisquared is the chi-squared merit function
% and q is the "goodness-of-fit" parameter.

% Fit the x y data to a straight line, ignoring 
% the uncertainties in x and y and determine the 
% uncertainty in the y direction due to the uncertainty
% in the x direction with the fitted slope.  This 
% works only if the data is very linear.
fit=polyfit(x,y,1);slope=fit(1);ysig2=slope*xsig;

% Combine the true uncertainty in y, ysig, with the 
% uncertainty in y determined from the x uncertainty, 
% ysig2, by using the square root of the sum of the squares.  
% This is only really vaild if the sources of error in 
% xsig and ysig are independent.
sigma=sqrt(ysig.^2 + ysig2.^2);sigma2=sigma.^2;

% Determine the best slope and y-intercept and 
% their uncertainties.  See Numerical Recipes, page 523.
S=sum(1./sigma2);
Sx=sum(xsig./sigma2);
Sx=sum(x./sigma2);
Sy=sum(y./sigma2);
Sxx=sum(x.*x./sigma2);
Sxy=sum(x.*y./sigma2);
Delta=S*Sxx-Sx*Sx;

a=(Sxx*Sy-Sx*Sxy)/Delta;     % a is the best fit y-intercept
b=(S*Sxy-Sx*Sy)/Delta;    % b is the best fit slope

asig=sqrt(Sxx/Delta);      % asig is the weighted uncertainty in a
bsig=sqrt(S/Delta);      % bsig is the weighted uncertainty in b

chisquared=sum(((y-a-b*x)./sigma).^2);    % the chi squared merit function
if chisquared ~= 0
  q=gammainc(length(x)/2-1,chisquared/2); % the goodness-of-fit parameter
else
  q=1;
end









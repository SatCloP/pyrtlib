% This function adds trailing blanks to fill a
% 1xn string corresponding to the number in input
% Usage:
%       [ss] = twodigstr(nn);
%       [ss] = twodigstr(nn,n);
% Input:
%        nn: number to convert in string
%         n: number of digit in output(default 2)
% Output:
%        ss: string of n digit corriponding to nn
% Example:
%        
%       ss = twodigstr(3);     % => ss = '03'
%       ss = twodigstr(30);    % => ss = '30'
%       ss = twodigstr(300,6); % => ss = '000300'
%       ss = twodigstr(3.2,6); % => ss = '0003.2'
%
% 2005/05/12 - Changed to include the option of number of digits

function [ss] = twodigstr(nn,n)

if nargin < 2
    n = 2;
end

ss = num2str(nn);
while length(ss) < n
   ss = ['0' ss];
end   

return
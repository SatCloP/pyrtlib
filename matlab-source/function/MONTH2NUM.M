% This function translate the name of a month in a number or viceversa.
%
% Input:
%        A number between 1 and 12 (1) or a string as 'JEN' (2)
%
% Output:
%        The string corresponding to numerical input (1) or the number corrisponding to input string (2)
% Es:
%    [n]=month2num('FEB'); or
%    [name]=month2num(12);

function [OUT]=month2num(INP)

months={'JAN';'FEB';'MAR';'APR';'MAY';'JUN';'JUL';'AUG';'SEP';'OCT';'NOV';'DEC'};

iclass = class(INP);

if iclass(1)=='d'
    
    OUT = char(months(INP));
    
elseif iclass(1)=='c'
    
    MON = upper(INP(1:3));
    LIS = strcmp(months,MON);
    OUT = find(LIS);
    
end    
    
return

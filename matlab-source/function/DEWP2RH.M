% This function converts dew point temperature (C) in relative humidity (RH).
%
% Inputs: 
%        Tdwp: Dew Point Temperature (C)
%        Tc:   Physical Temperature (C)
%
% Optional Input: 
%        ice:  If 1, then computes vapor pressure with respect to ice 
%              (default or 0 means with respect to water)
%
% Output: 
%        RH:  Relative Humidity (%)
%
% Es:
%    [RH] = Dewp2RH(Tdwp,Tc);


function [RH] = Dewp2RH(Tdwp,Tc,ice)

% for use in vapor pressure with respect to water
A =   7.5;
B = 237.3;

if nargin==3 & ice==1
   % for use in vapor pressure with respect to ice
   A =   9.5;
   B = 265.5;
end

% Vapor pressure (ED):
ED = 6.1078 * 10 .^ ((Tdwp * A)./(Tdwp + B)); % in mbs

% Saturated vapor pressure (ES):
ES = 6.1078 * 10 .^ ((Tc * A)./(Tc + B));

% Relative humidity = 100 * ED/ES
RH = 100 * ED./ES;

return


function gkg = ppmv2gkg(ppmv,gasid);

%
% function gkg = ppmv2gkg(ppmv,gasid);
%
% convert volume mixing ratio in ppmv to 
% mass mixing ratio in g/kg.
%
% inputs:
%   gkg:  mass mixing ratio (g/kg)
% gasid:  HITRAN gas id
%
% if gasid ==  1; H20 
% if gasid ==  2; CO2 
% if gasid ==  3; 03 
% if gasid ==  4; N2O 
% if gasid ==  5; CO 
% if gasid ==  6; CH4 
% if gasid ==  7; O2   
% if gasid ==  8; NO 
% if gasid ==  9; SO2 
% if gasid == 10; NO2 
% if gasid == 11; NH3 
% if gasid == 12; HNO3 
% if gasid == 13; OH 
% if gasid == 14; HF 
% if gasid == 15; HCL 
% if gasid == 16; HBR 
% if gasid == 17; HI 
% if gasid == 18; CLO 
% if gasid == 19; OCS 
% if gasid == 20; H2CO 
% if gasid == 21; HOCL 
% if gasid == 22; N2  
% if gasid == 23; HCN 
% if gasid == 24; CH3CL 
% if gasid == 25; H2O2 
% if gasid == 26; C2H2 
% if gasid == 27; C2H6 
% if gasid == 28; PH3 
% if gasid == 29; COF2 
% if gasid == 30; SF6 
% if gasid == 31; H2S 
% if gasid == 32; HCOOH 
% if gasid == 99; AIR 


% convert to parts per volume
ppv = ppmv / 1e6;

% multiply by ratio of masses to get mass 
% mixing ratio in g/g
gg = ppv * gas_mass(gasid)/gas_mass(99);

% multiply by 1000 to get g/kg
gkg = gg * 1000;


function mass_molecule = gas_mass(gasid);

% function mass_molecule = gas_mass(gasid);
%
% gas_mass.m returns the mass of the HITRAN gas ID gasid
%
% DCT 3/2/1996
% 
% NOTE: results are not accurate because amu 
% values need more significant figures
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


if gasid ==  1; amus = 2*1 + 16;  	end  	% H2O
if gasid ==  2; amus = 12 + 2*16;  	end   	% CO2
if gasid ==  3; amus = 3*16;		end 	% O3
if gasid ==  4; amus = 2*14 + 16;	end	% N2O
if gasid ==  5; amus = 12 + 16;  	end	% CO
if gasid ==  6; amus = 12 + 4*1;  	end	% CH4
if gasid ==  7; amus = 2*16;	  	end	% O2  
if gasid ==  8; amus = 14 + 16;  	end	% NO
if gasid ==  9; amus = 32 + 2*16;	end  	% SO2
if gasid == 10; amus = 14 + 2*16;  	end 	% NO2
if gasid == 11; amus = 14 + 3*1;  	end	% NH3
if gasid == 12; amus = 1 + 14 + 3*16;  	end	% HNO3
if gasid == 13; amus = 16 + 1;  	end	% OH
if gasid == 14; amus = 1 + 19;  	end	% HF
if gasid == 15; amus = 1 + 35;  	end	% HCL
if gasid == 16; amus = 1 + 80;  	end	% HBR
if gasid == 17; amus = 1 + 127;  	end	% HI
if gasid == 18; amus = 35 + 16;  	end	% CLO
if gasid == 19; amus = 16 + 12 + 32;  	end	% OCS
if gasid == 20; amus = 2*1 + 12 + 16;  	end	% H2CO
if gasid == 21; amus = 1 + 16 + 35;  	end	% HOCL
if gasid == 22; amus = 2*14;  		end	% N2 
if gasid == 23; amus = 1 + 12 + 14;  	end	% HCN
if gasid == 24; amus = 12 + 3*1 + 35;  	end	% CH3CL
if gasid == 25; amus = 2*1 + 2*16;  	end  	% H2O2
if gasid == 26; amus = 2*12 + 2*1;  	end	% C2H2
if gasid == 27; amus = 2*12 + 6*1;  	end	% C2H6
if gasid == 28; amus = 31 + 3*1;  	end	% PH3
if gasid == 29; amus = 12 + 16 + 2*19;  end 	% COF2
if gasid == 30; amus = 32 + 6*19;  	end	% SF6
if gasid == 31; amus = 2*1 + 32;  	end	% H2S
if gasid == 32; amus = 1 + 12 + 2*16 +1;end	% HCOOH
if gasid == 99; amus = 28.9402753668809;end % air (makes water/air=.621970585)

mass_proton = 1.6726485e-27;
mass_molecule = mass_proton*amus;



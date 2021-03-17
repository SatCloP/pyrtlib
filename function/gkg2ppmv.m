
function ppmv = gkg2ppmv(gkg,gasid);

%
% function ppmv = gkg2ppmv(gkg,gasid);
%
% convert mass mixing ratio in g/kg to 
% volume mixing ratio in ppmv.
%
% inputs:
%  g/kg:  mass mixing ratio (g/kg)
% gasid:  HITRAN gas id
%

% divide by 1000 to get g/g
gg = gkg / 1000;

% divide by ratio of masses to get volume
% mixing ratio in g/g
ppv = gg / (gas_mass(gasid)/gas_mass(99));

% convert to parts per million volume
ppmv = ppv * 1e6;

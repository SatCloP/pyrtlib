
function sonde = readsonde(filename);

%
% function sonde = readsonde(filename);
%
% load in a sonde and compute mixing ratios
% (currently set up to prompt for single sonde if filename not input)

% get filename if needed

if nargin == 0
  [filename,pathname] = uigetfile('/home/davet/WVIOP97/DATA/sgpsondewrpnC1.a1.*.cdf');
  filename = [pathname filename];
end

sonde = rd_netcdf(filename);

sonde.nss1970 = sonde.base_time+sonde.time_offset;
sonde.tdry = sonde.tdry+273.15;
sonde.dp = sonde.dp+273.15;
sonde.filename = filename;

% compute the mixing ratios from the P, T, and RHs
sonde.w1 = rh2mr(sonde.pres,sonde.tdry,sonde.rh);

% compute the mixing ratios from the P and DewPoints
sonde.w2 = dp2mr(sonde.pres,sonde.tdry,sonde.dp);


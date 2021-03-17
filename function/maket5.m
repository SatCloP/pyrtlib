
function maket5(wnbr1,wnbr2,pressure,temperature,gasID,JCHAR,ga)

%
% function maket5(wnbr1,wnbr2,pressure,temperature,gasID,JCHAR,ga)
%
% matlab function to produce homogenous path optical depth tape5 files.
% i.e. maket5(1400,1410,1013,296,2,'A',360)
%
% inputs:
% 	wnbr1,wnbr2  bounding wavenumbers (scalars)
%	pressure     total pressure in millibars (scalar)
%	temperature  temperature in K (scalar)
%	gasID	     HITRAN gasID of absorber species
%	JCHAR	     LBLRTM string character for absorber amount units
%	ga	     gas amount of absorber gasID in units of JCHAR
%

if ~exist('wnbr1') | isempty(wnbr1);wnbr1 = input('Enter minimum wavenumber : ');end
if ~exist('wnbr2') | isempty(wnbr2);wnbr2 = input('Enter maximum wavenumber : ');end
if ~exist('pressure') | isempty(pressure);pressure = input('Enter pressure [mbar] : ');end
if ~exist('temperature') | isempty(temperature);temperature = input('Enter temperature [K] : ');end
if ~exist('gasID') | isempty(gasID);gasID = input('Enter HITRAN gas ID : ');end
if ~exist('JCHAR') | isempty(JCHAR);JCHAR = input('Enter JCHAR, gas amount units flag : ','s');end
if ~exist('ga') | isempty(ga);ga = input('Enter absorber amount in JCHAR units : ');end

GA = zeros(7,1);
GA(gasID) = ga;

jchars = '     AA   CCCCCCC';
jchars(10+gasID) = JCHAR;

fid = fopen('tape5','w');

fprintf(fid,'%s\n','$ optical depth calculation for one homogenous path');
fprintf(fid,'%s\n',' HI=1 F4=1 CN=1 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AT=1 M=01 LS=0 MS=0 XS=0    0    0');
fprintf(fid,'%s\n','    0    1    1    1    1');
fprintf(fid,'%10.3f',wnbr1,wnbr2);fprintf(fid,'\n');
fprintf(fid,'%s\n','    0    1    1    1    1    7    1          0.000     0.000      .000    360.00');
fprintf(fid,'%s\n','     0.000                         0.001');
fprintf(fid,'%s\n','    1  one homogenous path');
fprintf(fid,'%10.3f',0,pressure,temperature);fprintf(fid,'%s\n',jchars);
fprintf(fid,'%10.3e',GA(1),GA(2),GA(3),GA(4),GA(5),GA(6),GA(7));fprintf(fid,'\n');
fprintf(fid,'%s\n','%');

fclose(fid);



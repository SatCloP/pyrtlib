
function writetape5(title,V1,V2,TBOUND,H1,H2,ANGLE,altitude,pressure,temperature,water_vapor,aflag)

%
% function writetape5(title,V1,V2,TBOUND,H1,H2,ANGLE,altitude,pressure,temperature,water_vapor,aflag)
%
% create a TAPE5 file
%
% this is specialized for uplooking calcs.  not bug free.
%
% things in {} correspond to variable names in the LBLRTM manual
%
% uses Gas IDs 1-7 and 4 X-sections
%
% user supplies altitude (km), pressure (mbar), temperature (K), and water vapor (g/kg)
% and the other gases (2-7) use aflag atmosphere.  
%
% (aflag =1, tropical; =2, midlatitude summer; =3, midlatitude winter; =4, subarctic summer;
%        =5, subarctic winter, =6, 1976 U. S. Standard)
%
% ---not cureently implimented--
% Automatically appends model atmosphere to 70 km w/ 2km thick layers. but cuts back number of layers
% to be == 200 if this exceeds 200 layers
% ---not cureently implimented--
% 
% Inputs: (all are required)
%
%       name    description				size       
% 1 	title	up to 80 characters of user info	up to 1x80
% 2	V1	min wavenumber bound			1x1
% 3	V2	max wavenumber bound			1x1
% 4	TBOUND	boundary surface (K)			1x1
% 5	H1	observer altitude (km)			1x1
% 6	H2	end point altitude (km)			1x1
% 7	ANGLE	zenith angle at H1 (degrees)		1x1
% 8	altitude layer boundary altitudes (km)		1xN
% 9	pressure layer boundary pressures (mbar)	1xN
% 10 	temperature layer boundary temperature (K)	1xN
% 11 	water_vapor layer boundary water vapor mixing ratio (g/kg)	1xN
% 12    aflag	flag for default atmosphere (1-6)
%
% Outputs:
%	 the file "tape5" is created
%

no_input_levels = length(altitude);

pressure = reshape(pressure,no_input_levels,1);
altitude = reshape(altitude,no_input_levels,1);
temperature = reshape(temperature,no_input_levels,1);
water_vapor = reshape(water_vapor,no_input_levels,1);

%% Append model atmosphere to input atmosphere to 70 km w/ 2km thick layers
%if max(altitude) <= 70
%	altitude = [altitude' fix(floor(max(altitude))/2)*2+2:2:70]';
%end

% keep number of levels <= 200
if length(altitude) > 200;altitude = altitude(1:200);end

no_total_levels = length(altitude);

fid = fopen('tape5','w');

% print out the title.  this line must begin with "$". {CXID}
if length(title) > 80;title = title(1:80);end
fprintf(fid,'$ %s \n',title);

% print out various parameters
%  HI {HIRAC}	=1 for Voigt profile
%  F4 {ILBLF4}	=1 lbl bound is 25 cm-1 for all layers
%  CN {ICNTNM}	=1 to include continua
%  AE {IAERSL}	=0 to exclude aerosols
%  EM {IEMIT}	=1 for radiance and transmittance calcs
%  SC {ISCAN}	=0 for no scanning function
%  FI {IFILTR}	=0 for no filter
%  PL {IPLT} 	=0 for no plotting
%  TS {ITEST}	=0 for no testing
%  AM {IATM}	=1 for LBLATM=yes
%  MG {IMRG}	=0 for normal merging
%  LA {ILAS}	=0 for no laser options
%  MS {IOD}	=0 for normal
%  XS {IXSECT}	=1 to include cross sections
%  {MPTS}	=0
%  {NPTS}	=0
fprintf(fid,' HI=1 F4=1 CN=1 AE=0 EM=1 SC=0 FI=0 PL=0 TS=0 AM=1 MG=0 LA=0 MS=0 XS=1   00   00\n',[]);

% print wavenumber bounds, V1 and V2 [cm^-1]
fprintf(fid,'%10.3e%10.3e\n',[V1 V2]);

% {TBOUND} boundary temperature
% print boundary parameters.  emissivity = 1;
fprintf(fid,'%7.3f   %5.3f    %5.3f    %5.3f\n',[TBOUND 1 0 0]);

% page 24 of manual for IATM=1
MODEL = 0; 	% user supplied atmospheric profile
ITYPE = 2; 	% slant path calculation
IBMAX = no_total_levels;	% number of user supplied layer boundaries
NOZERO = 1; 	% suppress zeroing absorber amounts
NOPRNT = 1; 	% selects short printout
NMOL   = 7; 	% number of molecules
IPUNCH = 0;
IFXTYP = 0;
MUNITS = 0;
RE = 0;
HSPACE = 100;
VBAR = (V1+V2)/2;
CO2MX = 360; 	% CO2 mixing ratio (ppmv)
fprintf(fid,'%5i%5i%5i%5i%5i%5i%5i%2i %2i%10.3f%10.3f%10.3f%10.3f\n',...
	[MODEL ITYPE IBMAX NOZERO NOPRNT NMOL IPUNCH IFXTYP MUNITS RE HSPACE VBAR CO2MX]);

% print slant path parameters
% {H1} observer altitude
% {H2} end point altitude
% {ANGLE} zenith angle at H1 =0 for uplooking = 180 for downlooking
fprintf(fid,'%10.3f%10.3f%10.3f\n',[H1 H2 ANGLE]);
for i = 1:no_total_levels
	fprintf(fid,'%10.3f',altitude(i));
	if rem(i,8)==0 & i ~= no_total_levels
		fprintf(fid,'\n',[]);
	end
end
fprintf(fid,'\n');

% print profile info
fprintf(fid,'%5i levels in the user defined profile\n',no_total_levels);
for i = 1:no_input_levels
 if aflag==1;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C111111\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==2;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C222222\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==3;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C333333\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==4;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C444444\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==5;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C555555\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==6;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C666666\n',[altitude(i) pressure(i) temperature(i)]);end
 fprintf(fid,'%10.3e%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',[water_vapor(i) 0 0 0 0 0 0]);
end
for i = no_input_levels+1:no_total_levels
 if aflag==1;fprintf(fid,'%10.3f%10.3f%10.3f     11   1111111\n',[altitude(i) 0 0]);end
 if aflag==2;fprintf(fid,'%10.3f%10.3f%10.3f     22   2222222\n',[altitude(i) 0 0]);end
 if aflag==3;fprintf(fid,'%10.3f%10.3f%10.3f     33   3333333\n',[altitude(i) 0 0]);end
 if aflag==4;fprintf(fid,'%10.3f%10.3f%10.3f     44   4444444\n',[altitude(i) 0 0]);end
 if aflag==5;fprintf(fid,'%10.3f%10.3f%10.3f     55   5555555\n',[altitude(i) 0 0]);end
 if aflag==6;fprintf(fid,'%10.3f%10.3f%10.3f     66   6666666\n',[altitude(i) 0 0]);end
 fprintf(fid,'%10.3e%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',[0 0 0 0 0 0 0]);
end

if 0
%
% add in the x-section info
%
fprintf(fid,'%5i%5i%5i selected x-sections are :\n',[13 0 0]);
fprintf(fid,'CCL4      CCL3F     CCL2F2    CCLF3     CHCL2F    CHCLF2    C2CL3F3\n');
fprintf(fid,'C2CL2F4   C2CLF5    CF4       CLONO2    HNO4      N2O5\n');

fprintf(fid,'%5i%5i  using U.S. Standard concentrations: \n',[no_total_levels 0]);
for i = 1:no_total_levels
	fprintf(fid,'%10.3f     6666666666666\n',altitude(i));
	fprintf(fid,'%10.3e',zeros(7,1));
	fprintf(fid,'\n');
	fprintf(fid,'%10.3e',zeros(6,1));
	fprintf(fid,'\n');
end
end

%
% add in the x-section info
%


fprintf(fid,'%5i%5i%5i selected x-sections are :\n',[4 0 0]);
fprintf(fid,'CCL4      CCL3F     CCL2F2    CHCLF2\n');

fprintf(fid,'%5i%5i  using U.S. Standard concentrations: \n',[no_total_levels 0]);
for i = 1:no_total_levels
	fprintf(fid,'%10.3f     6666\n',altitude(i));
	fprintf(fid,'%10.3e',zeros(4,1));
	fprintf(fid,'\n');
end

fprintf(fid,'%%\n');
fclose(fid);











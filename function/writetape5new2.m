
function writetape5new(title,VBOUND,TBOUND,HBOUND,zbound,altitude,pressure,temperature,water_vapor,aflag)

%
% function writetape5new(title,VBOUND,TBOUND,HBOUND,zbound,altitude,pressure,temperature,water_vapor,aflag)
%
% 

no_input_levels = length(altitude);
no_boundaries = length(zbound);

pressure = reshape(pressure,no_input_levels,1);
altitude = reshape(altitude,no_input_levels,1);
temperature = reshape(temperature,no_input_levels,1);
water_vapor = reshape(water_vapor,no_input_levels,1);

fid = fopen('tape5','w');

% title
fprintf(fid,'$ %s \n',title);

% flags
fprintf(fid,' HI=1 F4=1 CN=1 AE=0 EM=1 SC=0 FI=0 PL=0 TS=0 AM=1 MG=0 LA=0 MS=0 XS=1   00   00\n',[]);

% continuum multipliers
fprintf(fid,'%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',[1 1 1 1 1 1]);

% wavenumber bounds
fprintf(fid,'%10.3f%10.3f\n',[VBOUND(1) VBOUND(2)]);

% boundary temperature and emissivity
fprintf(fid,'%10.3f%10.3f\n',[TBOUND  1]);

MODEL = 0;			     % user supplied atmospheric profile
ITYPE = 2;			     % slant path calculation
IBMAX = no_boundaries;		     % number of user supplied layer boundaries
NOZERO = 1;			     % suppress zeroing absorber amounts
NOPRNT = 1;			     % selects short printout
NMOL   = 7;			     % number of molecules
IPUNCH = 0;
IFXTYP = 0;
MUNITS = 0;
RE = 0;
HSPACE = 100;
VBAR = (VBOUND(1)+VBOUND(2))/2;
CO2MX = 360; 	% CO2 mixing ratio (ppmv)
fprintf(fid,'%5i%5i%5i%5i%5i%5i%5i%2i %2i%10.3f%10.3f%10.3f%10.3f\n',...
	[MODEL ITYPE IBMAX NOZERO NOPRNT NMOL IPUNCH IFXTYP MUNITS RE HSPACE VBAR CO2MX]);

% print slant path parameters
% {H1} observer altitude = HBOUND(1)
% {H2} end point altitude = HBOUND(2)
% {ANGLE} zenith angle at H1 =0 for uplooking = 180 for downlooking = HBOUND(3)
fprintf(fid,'%10.3f%10.3f%10.3f\n',[HBOUND(1) HBOUND(2) HBOUND(3)]);
for i = 1:no_boundaries
	fprintf(fid,'%10.3f',zbound(i));
	if rem(i,8)==0 & i ~= no_boundaries
		fprintf(fid,'\n',[]);
	end
end
fprintf(fid,'\n');

% print profile info
fprintf(fid,'%5i levels in the user defined profile\n',no_input_levels);
for i = 1:no_input_levels
 if aflag==1;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C111111\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==2;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C222222\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==3;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C333333\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==4;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C444444\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==5;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C555555\n',[altitude(i) pressure(i) temperature(i)]);end
 if aflag==6;fprintf(fid,'%10.3f%10.3f%10.3f     AA   C666666\n',[altitude(i) pressure(i) temperature(i)]);end
 fprintf(fid,'%10.3f\n',[water_vapor(i)]);
end

%
% add in the x-section info
%

fprintf(fid,'%5i%5i%5i selected x-sections are :\n',[3 0 0]);
fprintf(fid,'CCL4      F11       F12 \n');
fprintf(fid,'%5i%5i  \n',[2 0]);
fprintf(fid,'%10.3f     AAA\n',min(altitude));
fprintf(fid,'%10.3e%10.3e%10.3e\n',[1.105e-04 2.783e-09 5.027e-04]);
fprintf(fid,'%10.3f     AAA\n',max(altitude));
fprintf(fid,'%10.3e%10.3e%10.3e\n',[1.105e-04 2.783e-09 5.027e-04]);

fprintf(fid,'%%\n');
fclose(fid);




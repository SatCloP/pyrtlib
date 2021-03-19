% This function runs AER LBLRTM with particular settings useful
% to compare simulations from RAOB with AERI (or FTIR) measurements.
%
% Usage:
%        [Rd,Tb,wn] = lblrtm(RAOB,NAME)
% Inputs:
%        RAOB: structure containing the following fields
%             YYYY  Year of launch
%             jtime Julian date-time
%             ZKm   Height [km]
%             Pmb   Pressure [mb]
%             TC    Temperature [C]
%             RH    Relative humidity [%]
%        NAME: string stating the range of wavenumbers. 
%              So far are available
%              'aerich1' NSA aeri channel 1 ([ 375.11 1799.85] cm-1)
%              'aerich2' NSA aeri channel 2 ([1799.85 3020.17] cm-1)
% Outputs:
%        Rd: computed radiance [ decaW/(m2 cm-1 sr) ]
%        Tb: computed brightness temperature [K]
%        wn: wavenumber [cm-1]
% Example:
%        load testraob
%        RAOB.TC = RAOB.TK - 273.16;
%        RAOB.ZKm= RAOB.ZKm- RAOB.ZKm(1);
%        RAOB.YYYY=1999; RAOB.jtime=60.5;
%        [Rd,Tb,wn] = lblrtm(RAOB,'aerich1');
%
% Nico, Sept, 2003
% Last Update 2003/10/09 - Version 8.1

function [Rd,Tb,wn] = lblrtm(RAOB,name)

%lblrundir = '/home/niko/univaq/LBLRTM/run/';      % Version 5.?
lblrundir = '/home/niko/SFTWR/LBLRTM/lblrtm/run/'; % Version 8.1

% purging existing tape files
if exist([lblrundir 'TAPE5']);   delete([lblrundir 'TAPE5']);   end;
if exist([lblrundir 'TAPE6']);   delete([lblrundir 'TAPE6']);   end;
if exist([lblrundir 'tape7']);   delete([lblrundir 'tape7']);   end;
if exist([lblrundir 'tape7tb']); delete([lblrundir 'tape7tb']); end;
if exist([lblrundir 'TAPE9']);   delete([lblrundir 'TAPE9']);   end;
if exist([lblrundir 'TAPE10']);  delete([lblrundir 'TAPE10']);  end;
if exist([lblrundir 'TAPE11']);  delete([lblrundir 'TAPE11']);  end;
if exist([lblrundir 'TAPE12']);  delete([lblrundir 'TAPE12']);  end;
if exist([lblrundir 'TMPX101']); delete([lblrundir 'TMPX101']); end;
if exist([lblrundir 'TMPX102']); delete([lblrundir 'TMPX102']); end;
if exist([lblrundir 'TMPX103']); delete([lblrundir 'TMPX103']); end;
if exist([lblrundir 'TMPX202']); delete([lblrundir 'TMPX202']); end;
if exist([lblrundir 'TMPX203']); delete([lblrundir 'TMPX203']); end;
if exist([lblrundir 'TMPXBIN']); delete([lblrundir 'TMPXBIN']); end;
if exist([lblrundir 'tp6']);     delete([lblrundir 'tp6']);     end;
if exist([lblrundir 'raobtmp.dat']); delete([lblrundir 'raobtmp.dat']); end;


% creating raob file
fid = fopen([lblrundir 'raobtmp.dat'],'w');
NL = length(RAOB.ZKm);
[MM,DD,hh,mm,ss] = jul2date(RAOB.YYYY,RAOB.jtime);
YYYYMMDD = [num2str(RAOB.YYYY) twodigstr(MM) twodigstr(DD)];
hhmmss = [twodigstr(hh) twodigstr(mm) twodigstr(round(ss))];
fprintf(fid,'%d %s %s\n',NL,YYYYMMDD,hhmmss);
for il = 1:NL
    fprintf(fid,'%8.4f %8.2f %8.2f %8.2f \n',RAOB.ZKm(il),RAOB.Pmb(il),RAOB.TC(il),RAOB.RH(il));
end
fclose(fid);

% moving to LBLRTM run directory
initdir = pwd;
cd(lblrundir);

% writing TAPE5, calling LBLRTM, converting output
% settape5_aeri is a special case made up for NSA1999 (Subarctic winter model. Look for H555555 in settape_aerich*.f)
eval(['! ./settape5_' name ' raobtmp.dat']);
! ./lblrtm
! ./prt11
! ./prt11tb

% retrieving output
MTX = load('tape7tb');

% clean up the mess
delete('TAPE5');
delete('TAPE6');
delete('tp6');
delete('tape7');
delete('tape7tb');
delete('TAPE9');
delete('TAPE10');
delete('TAPE11');
delete('TAPE12');
if exist('TMPX101'); delete('TMPX101'); end;
if exist('TMPX102'); delete('TMPX102'); end;
if exist('TMPX103'); delete('TMPX103'); end;
if exist('TMPX202'); delete('TMPX202'); end;
if exist('TMPX203'); delete('TMPX203'); end;
if exist('TMPXBIN'); delete('TMPXBIN'); end;
delete('raobtmp.dat');

% go back to initial directory
cd(initdir);

% outputs
wn = MTX(:,1);
Rd = MTX(:,2);
Tb = MTX(:,3);

return
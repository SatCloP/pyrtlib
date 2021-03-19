% Test for RTE Matlab functions wrt RTE fortran routines.
%
% Nico, Nov 2014
% Nico, Jan 2016 - Adapted for Mac
% Nico, Oct 2016 - Add clear/cloudy option

% Path settings
%path_Matools = 'C:\Matools\';
path_Matools = '/Users/Nico/Matools/';
addpath(path_Matools); addpath(mypath);
addpath([path_Matools 'RTE/']);
addpath([path_Matools 'RTMUli/']);

% Clear/Cloudy settings
cloudy = 0;

% Forward model absorption settings
%RTEFM = 'L87';         % ancora da convertire...
%RTEFM = 'L93';         % ancora da convertire...
RTEFM = 'rosen';        % differenze <=0.01 K (clear and cloudy)
%RTEFM = 'ros03';       % differenze <=0.01 K (clear and cloudy)
%RTEFM = 'Lilje04';     % ancora da convertire...
%MATFM.wvres = 'LIL04';
%MATFM.wvcnt = 'ROS98';
MATFM.wvres = 'ROS98'; % default
MATFM.wvcnt = 'ROS98';

% Loading
% fid = fopen([path_Matools 'RTE/test.dat']);
% nrb = 6;
% for i = 1:nrb
%   lin=[];  
%   while isempty(lin)
%     lin = fgetl(fid);
%   end
%   nlv = str2num(lin(10:11));
%   MTX = fscanf(fid,'%f',[4,nlv]);
%   MTX = MTX';  
%   info(i,1:40) = lin(1:40); 
%   RAOB(i).Pmb = MTX(:,1);
%   RAOB(i).ZKm = MTX(:,2);
%   RAOB(i).TK  = MTX(:,3) + 273.16;
%   RAOB(i).RH  = MTX(:,4);
%   clear MTX nlv; 
% end
% fclose(fid);
load([path_Matools 'RTE/test_raobs.mat']);
nrb = length(RAOB);

% Processing & printout
frq = [22.2400; 23.0400; 23.8400; 25.4400; 26.2400; 27.8400; 31.4000; 51.2600; 52.2800; 53.8600; 54.9400; 56.6600; 57.3000; 58.0000]';
%frq = [23.8, 31.6, 53.5, 55.5, 58.0];
ang = 90.;
nf = length(frq);

for i = 1:nrb
   
  z = RAOB(i).ZKm;
  p = RAOB(i).Pmb;
  t = RAOB(i).TK;
  rh= RAOB(i).RH;
  denliq = zeros(size(z)); % no cloud
  denice = zeros(size(z)); % no cloud
  cldh = zeros(2,0);
  if cloudy
     % build a cloud
     ib = 2; it = 4; denliq(ib:it) = 10*ones(it-ib+1,1); cldh(:,1) = [z(ib); z(it)]; 
     ib = 30; it = 32; denice(ib:it) = 0.1*ones(it-ib+1,1); cldh(:,2) = [z(ib); z(it)]; 
  end
  
  RTE = tbcloud(z,p,t,rh,denliq,denice,cldh,0,frq,ang,RTEFM);
  MAT = TbCloud_RTE(z,p,t,rh,denliq,denice,cldh,frq,ang,RTEFM);
  %MAT = TbCloud_Uli(z,p,t,rh,denliq,denice,cldh,frq,ang,MATFM);
  
  fprintf(1,['        PWV(cm) LWP(cm) IWP(cm)' repmat('%6.2f ',1,nf) ' \n'],frq);
  fprintf(1,['RTElib: %6.2f  %6.2f  %6.2f ' repmat('%6.2f ',1,nf) '\n'],RTE.srho,RTE.sliq,RTE.sice,RTE.tbtotal);
  fprintf(1,['RTEmat: %6.2f  %6.2f  %6.2f ' repmat('%6.2f ',1,nf) '\n'],MAT.srho,MAT.sliq,MAT.sice,MAT.tbtotal);
  fprintf(1,['Differ: %6.2f  %6.2f  %6.2f ' repmat('%6.2f ',1,nf) '\n'],RTE.srho-MAT.srho,RTE.sliq-MAT.sliq,RTE.sice-MAT.sice,RTE.tbtotal-MAT.tbtotal);
  
end  

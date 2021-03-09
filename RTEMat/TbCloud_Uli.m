% This uses routines adapted from Uli's idl code
%
% ... jcr ...........................................................
% User interface which computes brightness temperatures (Tb), mean
% radiating temperature (Tmr), and integrated absorption (Tau) for 
% clear or cloudy conditions,  Also returns all integrated quantities
% that the original TBMODEL, Cyber Version, returned.  The input
% profiles are not modified within this subroutine.  It is assumed
% that the input profiles start at the antenna height (zX(1)).  The
% input profiles must reach 50.0 mb.  This subroutine uses the
% algorithms described in Schroeder and Westwater (1991).
%
% INPUTS PASSED THRU ARGUMENT LIST: 
%   zX       - height profile (km MSL).
%   pX       - pressure profile (mb).
%   tkX      - temperature profile (K).
%   rhX      - relative humidity profile (fraction).
%   denliqX  - liquid density profile (g/m**3); density fraction = 1.0;
%   deniceX  - ice density profile (g/m**3); density fraction = 1.0;
%   frq      - channel frequencies (GHz).
%   anglesX  - elevation anglesX (deg).
%   nlayX    - number of cloud layers identified (<= maxc);
%   cldbaseX - cloud base heights (km MSL);
%   cldtopX  - cloud top heights  (km MSL);
%   iceX     - saturation vapor pressure switch.
%                iceX = 0: Compute saturation vapor pressure (es) 
%                          over water only.
%                iceX = 1: Compute es over water when temperature > -10 C.
%                          Compute es over ice when temperature <= -10 C.
%
%   OPTIONAL INPUTS
%   absmdl:      Absorption model for WV (default 'ROS98')
%                absmdl.wvres = wv resonant absorption
%                absmdl.wvcnt = wv continuuum absorption
%   Ray_tracing: if 1 (default) it computes ray tracing (RayTrac_xxx) for
%                distance between layers; otherwise use simple plane
%                parallel assumption (i.e. ds = diff(z)*airmass;
%
% OUTPUTS PASSED THRU ARGUMENT LIST:
%   tbtotalX - brightness temperature (K) includes cosmic background;
%              indexed by frequency and elevation angle
%   tbatmX   - atmospheric brightness temperature (K), no cosmic;
%              background;indexed by frequency and elevation angle
%   tmrX     - mean radiating temperature of the atmosphere (K);
%     	       indexed by frequency and elevation angle
%   tmrcldX  - mean radiating temperature (K) of the lowest cloud layer;
%              indexed by frequency and elevation angle
%   staudryX - dry air absorption integrated over each ray path (Np);
%              indexed by frequency and elevation angle
%   stauwetX - water vapor absorption integrated over each ray path (Np);
%              indexed by frequency and elevation angle
%   stauliqX - cloud liquid absorption integrated over each ray path (Np);
%              indexed by frequency and elevation angle
%   stauiceX - cloud ice absorption integrated over each ray path (Np);
%              indexed by frequency and elevation angle
%   srhoX    - water vapor density integrated along each ray path (cm);
%   	       indexed by elevation angle
%   swetX    - wet refractivity integrated along each ray path (cm);
%	       indexed by elevation angle
%   sdryX    - dry refractivity integrated along each ray path (cm);
%	       indexed by elevation angle
%   sliqX    - cloud ice density integrated along each ray path (cm);
%	       indexed by elevation angle
%   siceX    - cloud liquid density integrated along each ray path (cm);
%	       indexed by elevation angle
%
% RTE SUBROUTINES AND FUNCTIONS CALLED FROM HERE:
%   Bright_xxx   = compute temperature for the modified Planck radiance 
%   CldAbs_xxx   = computes cloud (liquid and ice) absorption profiles
%   CldInt_xxx   = integrates cloud water density of path ds (linear) 
%   CldTmr_xxx   = computes mean radiating temperature of a cloud 
%   ClrAbs_xxx   = computes clear-sky (h2o and o2) absorption profiles
%   ExpInt_xxx   = integrates (ln) absorption over profile layers
%   Planck_xxx   = computes modified planck radiance and related quantities
%   RayTrac_xxx  = computes refracted path length between profile levels
%   Refract_xxx  = computes vapor pressure and refractivity profiles
%   Vapor_xxx    = computes vapor pressure and vapor density 
% ...................................................................

%function [tbtotal,tbatm,tmr,tmrcld,tauwet,taudry,tauliq,tauice,srho,swet,sdry] = TbCloud_RTE(z,p,tk,rh,denliq,denice,frq,angles)
function OUT = TbCloud_RTE(z,p,tk,rh,denliq,denice,cldh,frq,angles,absmdl,Ray_tracing)

if nargin < 10
   absmdl.wvres = 'ROS98'; % default
   absmdl.wvcnt = 'ROS98'; % default
end

if nargin < 11
   Ray_tracing = 1; % default
end

% Settings
nl = length(z);
nf = length(frq);
nang = length(angles);
ncld = length(cldh(1,:));
ice = 0;

% Allocation
SPtaudry = zeros(nf,nang);
SPtauwet = zeros(nf,nang);
SPtauliq = zeros(nf,nang);
SPtauice = zeros(nf,nang);
Ptaudry = zeros(nf,nang,nl);
Ptauwet = zeros(nf,nang,nl);
Ptauliq = zeros(nf,nang,nl);
Ptauice = zeros(nf,nang,nl);
tbtotal = zeros(nf,nang);
tbatm = zeros(nf,nang);
tmr = zeros(nf,nang);
tmrcld = zeros(nf,nang);
srho = zeros(1,nang);
swet = zeros(1,nang);
sdry = zeros(1,nang);
sliq = zeros(1,nang);
sice = zeros(1,nang);

% ... convert height profile to (km above antenna height) ...
z0 = z(1);
z = z - z0;

% ... compute vapor pressure and vapor density ...
[e,rho] = Vapor_xxx(tk,rh,ice);

% ... convert cloud base and cloud top to (km above antenna height) ...
% ... compute (beglev) and (endlev) ...
cldh = cldh - z0;
for l = 1:ncld
    for i = 1:nl
        if z(i) == cldh(1,l); beglev(l) = i; end; 
        if z(i) == cldh(2,l); endlev(l) = i; end; 
    end
end

% ... compute refractivity ...
[dryn,wetn,refindx] = Refract_xxx(p,tk,e);

for k = 1:nang

    % ... Compute distance between each level (ds) ... 
    if Ray_tracing
       ds = RayTrac_xxx(z,refindx,angles(k),z0); % Ray tracing
    else 
       amass = 1 / sin( angles(k)*pi/180 ); % simple plane parallel
       ds = [0; diff(z)*amass];
    end
    % ds = [0; diff(z)]; % in alternative simple diff of z
    
    % ... Integrate over path (ds) ...
	srho(k) = ExpInt_xxx(1,rho,ds,1,nl,0.1);
	swet(k) = ExpInt_xxx(1,wetn,ds,1,nl,0.1);
	sdry(k) = ExpInt_xxx(1,dryn,ds,1,nl,0.1);
	if ncld > 0
	   sliq(k) = CldInt_xxx(denliq,ds,beglev,endlev);
	   sice(k) = CldInt_xxx(denice,ds,beglev,endlev);
    end

    
    % ... handle each frequency ...
    for j = 1:nf
        
        % ... compute the clear-air absorption profile ...
        [awet,adry] = ClrAbs_uli(p,tk,rho,frq(j),absmdl);  % here I use routines adapted from Uli's idl code

        % ... compute the cloudy-air absorption profile ...
        [aliq,aice] = CldAbs_xxx(tk,denliq,denice,frq(j));

        % ... Integrate awet over ds ... 
        [SPtauwet(j,k),Ptauwet(j,k,:)] = ExpInt_xxx(1,awet,ds,1,nl,1);
        % ... Integrate adry over ds ... 
        [SPtaudry(j,k),Ptaudry(j,k,:)] = ExpInt_xxx(1,adry,ds,1,nl,1);
        %  ... Integrate aliq over ds ... 
        [SPtauliq(j,k),Ptauliq(j,k,:)] = ExpInt_xxx(0,aliq,ds,1,nl,1);
        %  ... Integrate aice over ds ... 
        [SPtauice(j,k),Ptauice(j,k,:)] = ExpInt_xxx(0,aice,ds,1,nl,1);
        
        % ... Compute total layer-integrated absorption ...
        Ptaulay(j,k,:) = Ptauwet(j,k,:) + Ptaudry(j,k,:) + Ptauice(j,k,:) + Ptauliq(j,k,:);

        % ... Compute modified Planck radiance over related quantities ...
      	[boftotl,boftatm,boftmr,PSPtauprof,hvk] = Planck_xxx(frq(j),tk,Ptaulay(j,k,:));

        % ... Compute mean radiating temperature of lowest cloud layer ...
	    if ncld > 0
      	   tmrcld(j,k) = CldTmr_xxx(beglev(1),endlev(1),hvk,PSPtauprof,boftatm);
        end
        
        % ... assign output values ...
        tbtotal(j,k) = Bright_xxx(hvk,boftotl);
        tbatm(j,k) = Bright_xxx(hvk,boftatm(nl));
        tmr(j,k) = Bright_xxx(hvk,boftmr);
         
    end
   
end

OUT.tbtotal = tbtotal;
OUT.tbatm = tbatm;
OUT.tmr = tmr;
OUT.tmrcld = tmrcld;
OUT.tauwet = SPtauwet;
OUT.taudry = SPtaudry;
OUT.tauliq = SPtauliq;
OUT.tauice = SPtauice;
OUT.taulaywet = Ptauwet;
OUT.taulaydry = Ptaudry;
OUT.taulayliq = Ptauliq;
OUT.taulayice = Ptauice;
OUT.srho = srho;
OUT.swet = swet;
OUT.sdry = sdry;
OUT.sliq = sliq;
OUT.sice = sice;

return
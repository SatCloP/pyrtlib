% This was derived from NOAA RTE fortran routines
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

function OUT = TbCloud_RTE(z,p,tk,rh,denliq,denice,cldh,frq,angles,absmdl,Ray_tracing,Es)

if nargin < 10
   absmdl = 'rosen'; % default
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
    % this are based on NOAA RTE fortran routines
    for j = 1:nf

        switch absmdl
            case 'rosen' % default
                 % ... compute the clear-air absorption profile ...
                 [awet,adry] = ClrAbs_ros98(p,tk,e,frq(j));           
                 % ... compute the cloudy-air absorption profile ...
                 [aliq,aice] = CldAbs_ros98(tk,denliq,denice,frq(j));
                 
            case 'ros03'
                 % ... compute the clear-air absorption profile ...
                 [awet,adry] = ClrAbs_ros03(p,tk,e,frq(j));           
                 % ... compute the cloudy-air absorption profile ...
                 [aliq,aice] = CldAbs_ros03(tk,denliq,denice,frq(j));
                                
            case 'ros16'
                % using routines provided by P. Rosenkranz on Aug 10 2016
                % this is the same as Mak11 (I keep it for hystorical reasons)
                [awet,adry] = ClrAbs_ros16(p,tk,e,frq(j));   
                [aliq,aice] = CldAbs_ros16(tk,denliq,denice,frq(j));

            case 'Mak11'
                % using routines provided by P. Rosenkranz on Aug 10 2016
                % (using o2abs_d7.f, i.e. Makarov et al., 2011)
                [awet,adry] = ClrAbs_Mak11(p,tk,e,frq(j));          
                [aliq,aice] = CldAbs_ros16(tk,denliq,denice,frq(j));

            case 'ros17'
                % Rosenkranz, P.W.: Line-by-line microwave radiative transfer (non-scattering), Remote Sens. Code Library, doi:10.21982/M81013, 2017
                [awet,adry] = ClrAbs_ros17(p,tk,e,frq(j));   
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16
                
            case 'ros18'
                % Rosenkranz, personal communication, 2018/06/20 (email)
                [awet,adry] = ClrAbs_ros18(p,tk,e,frq(j));   
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16

            case 'ros18sd'
                % Rosenkranz, personal communication, 2018/12/13 (email)
                [awet,adry] = ClrAbs_ros18_SD(p,tk,e,frq(j)); % only at 22GHz! To be implemented for 118 GHz from o2abs_1-sd.f
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16

            case 'ros19'
                % Rosenkranz, personal communication, 2019/03/18 (email)
                [awet,adry] = ClrAbs_ros19(p,tk,e,frq(j)); % New coefficients but no SD! 
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16
                                
            case 'ros19sd'
                % Rosenkranz, personal communication, 2019/02/12 (email)
                [awet,adry] = ClrAbs_ros19_SD(p,tk,e,frq(j)); % only at 22 and 183GHz! To be implemented for 118 GHz from o2abs_1-sd.f
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16

            case 'ros20'
                % Rosenkranz, personal communication, 2020/02/27 (email)
                [awet,adry] = ClrAbs_ros20(p,tk,e,frq(j)); % Updated O2 abs wrt 19 
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16
                                
            case 'ros20sd'
                % Rosenkranz, personal communication, 2020/02/27 (email)
                [awet,adry] = ClrAbs_ros20_SD(p,tk,e,frq(j)); % Updated O2 abs wrt 19 (to be implemented for 118 GHz from o2abs_1-sd.f)
                [aliq,aice] = CldAbs_ros17(tk,denliq,denice,frq(j)); % this is the same as ros16
                                
            case 'MPM18'
                % da implementare starting from ros17 
                %(considerando Tretyakov et al. 2005; Makarov et al., 2011; Makarov et al., 2013; Koshelev et al., 2015; Koshelev et al., 2016; Tretyakov 2016)
                
            case 'uncertainty'
                [awet,adry] = ClrAbs_uncertainty(p,tk,e,frq(j));
                % da implementare (considerando Rosenkranz 2015, vedi Phil's routine abliq12.f)
                %[aliq,aice] = CldAbs_uncertainty(tk,denliq,denice,frq(j));
                [aliq,aice] = CldAbs_ros16(tk,denliq,denice,frq(j)); % da sostituire con riga sopra!

            case 'uncertainty_r18'
                [awet,adry] = ClrAbs_uncertainty_r18(p,tk,e,frq(j));
                % da implementare (considerando Rosenkranz 2015, vedi Phil's routine abliq12.f)
                %[aliq,aice] = CldAbs_uncertainty(tk,denliq,denice,frq(j));
                [aliq,aice] = CldAbs_ros16(tk,denliq,denice,frq(j)); % da sostituire con riga sopra!
                
            case 'uncertainty_r20sd'
                [awet,adry] = ClrAbs_uncertainty_ros20sd(p,tk,e,frq(j)); % Updated O2 abs wrt 19 (to be implemented for 118 GHz from o2abs_1-sd.f)
                % da implementare (considerando Rosenkranz 2015, vedi Phil's routine abliq12.f)
                %[aliq,aice] = CldAbs_uncertainty(tk,denliq,denice,frq(j));
                [aliq,aice] = CldAbs_ros16(tk,denliq,denice,frq(j)); % da sostituire con riga sopra!
                                
            case 'mix' % here one can mix models up (e.g. Ros98 for clear, Ros16 for liq)
                [awet,adry] = ClrAbs_ros98(p,tk,e,frq(j));           
                [aliq,aice] = CldAbs_ros16(tk,denliq,denice,frq(j));
               
        end
            

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
      	%[boftotl,boftatm,boftmr,PSPtauprof,hvk] = Planck_xxx(frq(j),tk,Ptaulay(j,k,:));
      	[boftotl,boftatm,boftmr,PSPtauprof,hvk] = Planck_SAT(frq(j),tk,Ptaulay(j,k,:),Es);

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
% JOE SHAW ROUTINE
% Water_index computes infrared complex refractive 
% index, reflectance and emittance of water, using 
% the refractive index from "Optical constants of 
% water in the 200-nm to 200-micron wavelength region," 
% G.M. Hale and M.R. Querry, Appl. Opt. 12(3), 555-563, 
% March 1973.
% The Fresnel calculation used here explicitly 
% recognizes the dependence on both wavelength and 
% angle so that the correct answer results when
% multiple wavelengths are used.
%
% N.B.: L'ho testato; per ang=45 e wvlng=14.2 da gli stessi risultati di dielectriconstant

wl = (3:0.1:15)';
angles = [45];
num_angles=length(angles);
load wat_nk_HQ.txt;
wl_coarse = wat_nk_HQ(:,1);
n_coarse = wat_nk_HQ(:,3);
k_coarse = wat_nk_HQ(:,4);
dn_coarse = wat_nk_HQ(:,5);
dk_coarse = wat_nk_HQ(:,6);
clear wat_nk_HQ;

% ...interpolate wavelength scale...
n = spline(wl_coarse,n_coarse,wl);
k = spline(wl_coarse,k_coarse,wl);
n_salt = spline(wl_coarse,(n_coarse+dn_coarse),wl);
k_salt = spline(wl_coarse,(k_coarse+dk_coarse),wl);
clear wn_coarse n_coarse k_coarse dn_coarse dk_coarse;
nk = n-k*i;
plot(wl,n,wl,k+1);
pause;
%nk_salt = n_salt-k_salt*i;
%nk = n_salt-k_salt*i;

%...reflectances...
theta = angles.*(2.0*pi/360.0);
for m=1:num_angles;
   th_ref(:,m) = asin(sin(theta(m))./nk);
end;

for m=1:num_angles;
   Rp(:,m) = (nk.*cos(theta(m))-cos(th_ref(:,m))) ./ (nk.*cos(theta(m))+cos(th_ref(:,m)));
   Rs(:,m) = (cos(theta(m))-nk.*cos(th_ref(:,m))) ./ (cos(theta(m))+nk.*cos(th_ref(:,m)));
   RRp(:,m) = Rp(:,m).*conj(Rp(:,m));
   RRs(:,m) = Rs(:,m).*conj(Rs(:,m));
end;

plot(wl,RRp, wl,RRs);
xlabel('wavelength (microns)'),...
%plot(angles,RRp, angles,RRs);
%xlabel('angle (deg)'),...
ylabel('Reflectance'),...
title('p and s reflectances');
pause;

%...emittances...
for m=1:num_angles;
   ep(:,m) = 1-RRp(:,m);
   es(:,m) = 1-RRs(:,m);
end;

plot(wl,ep, wl,es);
xlabel('wavelength (microns)'),...
%plot(angles,ep, angles,es),...;
%xlabel('angle (deg)'),...   
ylabel('Emittance'),...
title('p and s emittances');
pause;

%clear m k n th_ref Rp Rs;
% Function dielectriconstant
%
% This function calculates the complex dielectric constant of sea water 
% in the microwave or in a infrared region (13.6-14.8 micron) of the spectra, 
% giving in input salinity (psu), sea temperature (C), wavelength (cm).
%
% Es:
%    [eps]=dielectriconstant(slnty,temp,wvlength,flag)
%
% INPUTS:
%         slnty: water salinity (psu (?))     (unused if flag='ir')
%         temp: sea water temperature         (unused if flag='ir')
%         wvlength: wave length in cm
%         flag: 'mw' for microwave, 'ir' for infrared
%
% NOTES:
%        The function [eps,rH,rV]=epsilon(sol,temp,wave,tht) 
%        in C:/Nico/5mm/mat5mm is now obsolete
%
% REFERENCES:
%
%  Microwave: 
%             Stogryn, M.
%             "Equations for calculating the dielectic constant of saline water"
%             IEEE Trans. Trans. Microwave Theory and Tech., V19, N8, 733-736, 1971
%           
%   Infrared: G.M. Hale and M.R. Querry,
%             "Optical constants of water in the 200-nm to 200-micron wavelength region," 
%             Appl. Opt. 12(3), 555-563, March 1973.
%
% Nico, 8/2000

function [eps,n,k]=dielectriconstant(slnty,temp,wvlength,flag)

  switch flag

   case 'mw'
      
%    Stogryn, M.
%    "Equations for calculating the dielectic constant of saline water"
%    IEEE Trans. Trans. Microwave Theory and Tech., V19, N8, 733-736, 1971

      nopm = slnty*0.01*(1.707+0.001*(1.205*slnty+0.0004058*slnty^2));
      sig = nopm*(10.394-2.3776*nopm+0.68258*nopm^2-0.13538*nopm^3+0.010086*nopm^4);
      d = 25. - temp;
      cond = sig*(1-0.01*(1.962*d-0.00808*d^2)-d*nopm*(0.00001*(3.02+3.922*d)+nopm*0.00001*(1.721-0.6584*d)))*0.01;
      ls = 3.*(1.1109-0.03824*temp+0.0001*(6.938*temp^2-0.05096*temp^3));
      es = 87.74-0.4008*temp+0.0001*(9.398*temp^2+0.0141*temp^3);
      eh = 4.9;
      zn = 1.+(ls/wvlength)^2;
      epsr = eh+(es-eh)/zn;
      epsi = (es-eh)*(ls/wvlength)/zn;
      epsi = epsi+60.*wvlength*cond;
      eps=epsr+j*epsi;
      n=epsr;
      k=epsi;
      
   case 'ir'
      
      % N.B.: OLD version; works just for wavelegth=14 microns
      
      %w1 = wvlength*1.e4;;
      %w2 = w1^2;
      %w3 = w1^3;    
      %n_r = -0.005148965*w3 + 0.2141723*w2 - 2.9029145*w1 + 14.001469339;
      %n_i = 0.001834057*w3 - 0.09327748*w2 + 1.5790118*w1 - 8.4867982;
      %cindex = n_r - j * n_i;
      %eps = cindex^2;
      
      % N.B.: NEW version; works for wavelegth from 3 up to 15 microns
      
      % It computes infrared complex refractive index using the 
      % refractive index from "Optical constants of water in the 
      % 200-nm to 200-micron wavelength region," G.M. Hale and 
      % M.R. Querry, Appl. Opt. 12(3), 555-563, March 1973.
      
      load wat_nk_HQ.txt;
      wl_coarse = wat_nk_HQ(:,1);
      n_coarse = wat_nk_HQ(:,3);
      k_coarse = wat_nk_HQ(:,4);
      % ...interpolate wavelength scale...
      wl = wvlength*1.e4; % cm -> microns
      n = spline(wl_coarse,n_coarse,wl);
      k = spline(wl_coarse,k_coarse,wl);
      nk = n - k*i;
      eps = nk^2;
      
  end


return
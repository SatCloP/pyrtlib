% This routine will provide values and units for all the 
% universal constants that I needed in my work.
%
% Reference (not for all): 
% P.J. Mohr, B.N. Taylor, and D.B. Newell (2015), "The 2014 CODATA Recommended
% Values of the Fundamental Physical Constants" (Web Version 7.0), http://physics.nist.gov/cuu/index.html
% Values as of 11/12/2015
% 
% Usage:
%       [const,units]=constants(string);
% Input:
%       string: String specifying which constant is needed
%               Avalaible right now:
%               'avogadro'            Avogadro number [mol-1] 
%               'boltzmann'         Boltzmann constant [J K-1]  
%               'EarthRadius'       Earth radius [km]
%               'light'             Light speed [m s-1]  
%               'Np2dB'             Neper to Decibel [dB/Np]  
%               'planck'            Planck constant [J Hz-1] 
%               'Rdry'              Gas constant of dry air [J kg-1 K-1]
%               'Rwatvap'           Gas constant of water vapor [J kg-1 K-1]
%               'Tcosmicbkg'        Cosmic Background Temperature [K]
% Outputs:
%         const: Numerical Value of the asked constant
%         units: String specifying which units are used
% Es:   
%       [h,units]=constants('planck');
%
% Nico, 2000
% 
% History
% 2003/06/06 - Added 'Rdry' and 'Rwatvap'
% 2005/12/10 - Added 'Tcosmicbkg'
% 2021/03/09 - Added 'avogadro', 'gravity'

function [out,units]=constants(strng)

switch strng
   case 'avogadro'
      nA = 6.022140857E23; % [mol-1] Avogadro's number
      units='[mol-1]';
      out= nA;
   case 'boltzmann'
      K=1.380658 * 1e-23; % [J K-1] Boltzmann constant
      units='[J K-1]';
      out=K;
   case 'EarthRadius'
      R=6370.949; % [Km] Earth radius
      units='[Km]';
      out=R;
   case 'gravity'
      g = 9.80665; % [m s-2] Acceleration due to gravity 
      units='[m s-2]';
      out=g;
   case 'light'
      c=299792458; % [m s-1] Light speed
      units='[m s-1]';
      out=c;
   case 'Np2dB'
      Np2dB=10*log10(exp(1));
      units='[dB/Np]';
      out=Np2dB;
   case 'planck'
      h=6.6260755 * 1e-34; % [J Hz-1] Planck constant 
      units='[J Hz-1]';
      out=h;
   case 'Rdry'
      Rd=287.04; % [J kg-1 K-1] Gas constant of dry air
      units='[J kg-1 K-1]';
      out=Rd;
   case 'Rwatvap'
      Rv=461.50; % [J kg-1 K-1] Gas constant of water vapor 
      units='[J kg-1 K-1]';
      out=Rv;
   case 'Tcosmicbkg'
      %Tcos = 2.736; % +/- 0.017 [K] Cosmic Background Temperature, from Janssen, Atmospheric Remote Sensing by Microwave Radiometry, pag.12
      Tcos = 2.728; % +/- 0.004 [K] Cosmic Background Temperature, from Fixsen D.J. et al 1996 (Astrophysical Journal v.473, p.576), from full COBE FIRAS Data Set
      units='[K]';
      out=Tcos;
     
   otherwise
      fprintf(1,' No constant avalaible with this name: "%s" . Sorry... \n',strng)
      out=[];
      units=[];
      
end
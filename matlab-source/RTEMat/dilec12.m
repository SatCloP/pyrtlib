% c   Purpose: Computes the complex dielectric constant for liquid water,
% c   with a negative imaginary part representing dissipation.
% 
% c   Complex logarithm is used here. It should be defined with
% c   imaginary part in the range -pi to +pi.
% c
% c   Copyright ? P.W. Rosenkranz  Apr. 15, 2014
% c   Creative Commons license CC BY-SA
%
% kappa = dilec12(f,tk)
%
%c     inputs:
%      real f  ! frequency in GHz, 
%      real tk ! Kelvin temperature
%c     validated for 20<f<220 GHz at 248<tk<273; 1<f<1000 GHz at 273<tk<330.
%
%c     outputs:
%      complex kappa ! complex dielectric constant
%
% References:
%
%c  static dielectric constant model from
%c  Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
%c  Debye term from 
%c  W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
%c  B band from
%c  P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. v.53(3) pp.1387-93 (2015).

function kappa = dilec12(f,tk)

%      implicit none
%c   arguments-

%c   local variables: 
%      real sd ! Debye freq.
%      real tc,hdelta,delta,f1,theta
%      complex chip,chij,dchi,cnorm,z,z1,z2

      tc = tk - 273.15;
      z = complex(0.,f);
      theta = 300./tk;

%c  static dielectric constant model from
%c  Patek et al. (J.Phys.Chem.Ref.Data. v.38(1), 21 (2009).
      kappa = -43.7527*theta^.05 + 299.504*theta^1.47 ... 
             -399.364*theta^2.11 + 221.327*theta^2.31;

%c  Debye term from 
%c  W. Ellison, J. Phys. Chem. Ref. Data, 36, 1-18 (2007).
      delta = 80.69715 * exp(-tc/226.45);
      sd = 1164.023 * exp( -651.4728 / (tc+133.07) );
      kappa = kappa - delta*z / (sd+z);

%c  B band from
%c   P.W. Rosenkranz, IEEE Trans. Geosci. & Remote Sens. 
%c   v.53(3) pp.1387-93 (2015).
      delta = 4.008724 * exp(-tc/103.05);
      hdelta = delta/2.;
      f1 = 10.46012 + 0.1454962*tc + 6.3267156E-02 * tc^2 + 9.3786645E-04*tc^3;
%       z1 = (-.75,1.) * f1;
%       z2 = (-4500.,2000.)
      z1 = complex(-.75,1.) * f1;
      z2 = complex(-4500.,2000.);
      cnorm = log(z2/z1);
      chip = hdelta * log( (z-z2)/(z-z1) ) / cnorm;
      chij = hdelta * log( (z-conj(z2))/(z-conj(z1)) ) / conj(cnorm);
      dchi = chip + chij - delta;
      kappa = kappa + dchi;

      return
      end


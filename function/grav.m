%=======================================================================
%
% function GRAV = grav(Z, WINDE, WINDN, LAT, LON)
%
%=======================================================================
%
%              University of Maryland Baltimore County [UMBC]
%
%              AIRS
%
%              GRAV
%
%
% ROUTINE NAME: GRAV
%
% ABSTRACT: Function for computing Earth's gravity.
%
% CALL PROTOCOL:
%    GRAV(Z, WINDE, WINDN, LAT, LON)
%
% INPUT PARAMETERS:
%    type      name    purpose                     units
%    --------  ------  --------------------------  ---------------------
%    REAL      LAT     latitude                    degrees          
%    REAL      LON     longitude                   degrees
%    REAL      WINDE   wind velecity east          m/s
%    REAL      WINDN   wind velocity north         m/s
%    REAL      Z       altitude                    m
%
% OUTPUT PARAMETERS:
%    type      name    purpose                     units
%    --------  ------  --------------------------  ---------------------
%    REAL      GRAV    Earth gravity               m/s^2
%
%    Function to calculate Earth gravity (gravitation plus
%    centripetal acceleration) for points in the atmosphere.
%
%    It calculates surface gravity using an equation given in
%    "American Institute of Physics Handbook", 1963, page 2-102.
%    This equation is essentially a variant of the International
%    Gravity Formula with an extra term for longitude (which is
%    very minor).
%
%    Centripetal acceleration is tangental velocity squared over the
%    radius.
%
%    Gravitation is the gravitational constant times the Earth's mass
%    divided by the square of the radius.
%
%    Gravity at any point in the atmosphere is simply surface gravity
%    plus the change in gravitation and centripetal acceleration at
%    that point compared to the surface.
%

function GRAV = grav(Z, WINDE, WINDN, LAT, LON);

%      =================================================================
%      REAL FUNCTION GRAV(Z, WINDE, WINDN, LAT, LON)
%      =================================================================

%      Constants for (1 - b^2/a^2) with
%      a = 6.378388E+6 m = equatorial radius, and
%      b = 6.356911E+6 m = polar radius.
%
%      Constants for normal gravity equation
%      (see "American Institute of Physics Handbook", 1963, pg 2-102) 
%
%      Constants for pi/180, 2*pi, and Earth's rotational speed
%      of w=1/86400 rev/s

B2 = 4.041031E+13;
ABTERM = 6.724285E-3;
G0 = 9.780455;
C1 =5.30157E-3;
C2 = -5.85E-6;
C3 = 6.40E-6;
PI180 = 1.7453293E-2;
PI2 = 6.28318531;
W = 1.1574074E-5;

%    Calculate longitude term
%    Add offset of 18 degrees, double it, convert to radians, and
%    take the cosine
       COSLON = cos(PI180*2.0*(LON + 18.0));

%    Calculate the latitude terms
%    Convert Latitude into radians
       LTRAD=PI180*LAT;

%    Calculate sine and cosine terms
       COSLT = cos(LTRAD);
       COSLT2 = COSLT.^2;
       SINLT2 = ( sin(LTRAD ) ).^2;
       SIN2LT = ( sin( 2.0*LTRAD ) ).^2;

%    Calculate the Earth's radius at this latitude
       R = sqrt( B2/( 1.0 - COSLT2*ABTERM ) );

%    Calculate total distance from Earth's center
       RTOT = R + Z;

%    Calculate gravity at the Earth's surface
       G_SUR = G0*( 1.0 + C1*SINLT2 + C2*SIN2LT + C3*COSLT2*COSLON );

%    Calculate the centripetal term at the Earth's surface
       C_SUR = COSLT2*R*(PI2*W).^2;

%    Calculate the centripetal term at altitude z (with wind)
       C_Z = ( ( PI2*RTOT*COSLT*W + WINDE ).^2 + (WINDN).^2 )/RTOT;

%    Calculate the change in gravitation with altitude
       GRAVZ=(G_SUR + C_SUR)*(1.0 - R.^2/RTOT.^2);

       GRAV = G_SUR + (C_SUR - C_Z) - GRAVZ;

% ...................................................................
%     Ray-tracing algorithm of Dutton, Thayer, and Westwater, rewritten for
%     readability & attempted documentation.  Based on the technique shown in
%     Radio Meteorology by Bean and Dutton (Fig. 3.20 and surrounding text).
%
%     inputs passed as arguments: 
%          z       = height profile (km above observation height, z0) 
%          refindx = refractive index profile
%          angle   = elevation angle (degrees)
%          z0      = observation height (km msl)
%     outputs: 
%          ds      = array containing slant path length profiles (km) 
% ...................................................................

function ds = RayTrac_xxx (z,refindx,angle,z0);

deg2rad = pi / 180;
re = constants('EarthRadius'); % Earth radius [km]
nl = length(z);

% Check for refractive index values that will blow up calculations.
for i = 1:nl
    if refindx(i) < 1
       display('RayTrac_xxx: Negative rafractive index');
       return
    end
end

% If angle is close to 90 degrees, make ds a height difference profile.
if (angle >= 89 & angle <= 91) | (angle >= -91 & angle <= -89) 
    ds(1) = 0.0;
    for i = 2:nl
        ds(i) = z(i) - z(i-1);
    end
end

% The rest of the subroutine applies only to angle other than 90 degrees.
% Convert angle degrees to radians.  Initialize constant values.

theta0 = angle * deg2rad;
rs = re  +  z(1) +  z0;
costh0 = cos(theta0);
sina = sin(theta0 * 0.5);
a0 = 2. * (sina^2);
  
% Initialize lower boundary values for 1st layer.
ds(1)= 0.;
phil = 0.;
taul = 0.;
rl = re  +  z(1) +  z0;
tanthl = tan(theta0);


% Construct the slant path length profile.
for i = 2:nl
    
    r = re + z(i) + z0;
    
    % Compute layer-average refractive index.
    if ( refindx(i)==refindx(i-1) | refindx(i)==1. | refindx(i-1) == 1. )
        refbar = (refindx(i) + refindx(i-1)) * 0.5;
    else
        refbar = 1. + (refindx(i-1) - refindx(i)) / (log ((refindx(i-1) - 1.) / (refindx(i) - 1.)));
    end
        
    argdth = z(i) / rs - ( (refindx(1) - refindx(i) ) * costh0 / refindx(i) );
    argth = 0.5 * (a0 + argdth) / r;
        
    % Check for ducting. If found, set invalid ray path flag + skip out.
    if argth <= 0 
       display(['RayTrac_xxx: Ducting at ' num2str(angle) 'degrees']);
       return
    end
            
    % Compute d-theta for this layer.
    sint = sqrt (r * argth);
    theta = 2. * asin (sint);
    if (theta - 2. * theta0) <= 0.
        dendth = 2. * (sint + sina) * cos((theta + theta0) * 0.25);
        sind4 = (0.5 * argdth - z(i) * argth) / dendth;
        dtheta = 4. * asin (sind4);
        theta = theta0 + dtheta;
    else
        dtheta = theta - theta0;
    end
                
    % Compute d-tau for this layer (eq.3.71) and add to integral, tau.
    tanth = tan (theta);
    cthbar = ((1. / tanth) + (1. / tanthl)) * 0.5;
    dtau = cthbar * (refindx(i-1) - refindx(i)) / refbar;
    tau = taul + dtau;
    phi = dtheta + tau;
    
    % Compute length of arc across this layer (segment qf in Fig.3.20).
    ds(i) = sqrt( (z(i) - z(i-1))^2 + 4. * r * rl * ((sin((phi - phil) * 0.5))^2));
    if dtau ~= 0.
       dtaua = abs(tau - taul);
       ds(i) = ds(i) * (dtaua / (2.* sin (dtaua * 0.5)));
    end
    
    % Make upper boundary into lower boundary for next layer.
    phil = phi;
    taul = tau;
    rl = r;
    tanthl = tanth;
               
end
            


return

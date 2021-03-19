% Airmass computes air mass giving elevation angles (degrees) or inverse.
% Giving in input the absorption scale height (Km) it accounts for earth 
% curvature, following Han and Westwater, 2000 (Eq. 28).
% Giving a string, it uses the model in Meunier et al, 2013 (Eq.6-7)
% 
% Es:
%    [amass]=airmass(el_angles);     % secant model
% or 
%    [amass]=airmass(el_angles,1.8); % model in Han and Westwater, 2000 (Eq. 28)
% or 
%    [amass]=airmass(el_angles,'M'); % model in Meunier et al, 2013 (Eq.6-7)
% or 
%    [el_angles]=airmass(amass,1,'r'); % secant model
%
% Nico,2000
% Nico,2019 - modifies to add the model in Meunier et al, 2013 (Eq.6-7)

function [out]=airmass(varargin)

if nargin == 1
   
   ang = varargin{1};
   amass = 1 ./ sin(pi*ang/180.);  
   out=amass;
   
elseif nargin == 2
   
   ang = varargin{1};
   H = varargin{2}; % Water vapor scale height [Km]
   R = constants('EarthRadius'); % Earth radius [Km]
   
   if isnumeric(H)
      % Han and Westwater, 2000 (Eq. 28)  
      amass0 = 1 ./ sin(pi*ang/180.);     
      amass = amass0 - H * amass0 .* (amass0.^2 - 1.)/R;
      out=amass;
   else
      % Meunier et al, 2013 (Eq.6-7)
      ae = R * 4/3;
      z_km = glatm(1);
      nlev = length(z_km);
      nang = length(ang);
      out = zeros(size(ang));
      for ia = 1:nang
          ltot = 0;
          phi = ang(ia) * pi/180;
          for il = 1:nlev-1
              dz = z_km(il+1)-z_km(il);
              aesinphi = ae * sin(phi);
              l = sqrt(aesinphi^2 + dz^2 + 2*ae*dz) - aesinphi;
              phi = asin( l*cos(phi)/ae + sin(phi));
              ltot = ltot + l;
           end
           out(ia) = ltot/z_km(end);
      end
   end
   
elseif nargin == 3  % Inverse: from amass to ang
   
   amass = varargin{1};
   ang = asin(1./amass);
   ang = ang * 180/pi ;   
   out=ang;
   
end

return
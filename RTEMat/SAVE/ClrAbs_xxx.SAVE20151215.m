%   Computes profiles of water vapor and dry air absorption for
%   a given set of frequencies.  Subroutines H2O_xxx and O2_xxx
%   contain the absorption model of Leibe and Layton [1987: 
%   Millimeter-wave properties of the atmosphere: laboratory studies 
%   and propagation modeling. NTIA Report 87-224, 74pp.] with oxygen 
%   interference coefficients from Rosenkranz [1988: Interference 
%   coefficients for overlapping oxygen lines in air.  
%   J. Quant. Spectrosc. Radiat. Transfer, 39, 287-97.]
%
%     inputs passed as arguments: 
%          p      = pressure profile (mb)
%          tk     = temperature profile (K)
%          !!!e      = vapor pressure profile (mb)
%          rho    = water vapor density (g/m3)
%          frq    = frequency (GHz)
%       absmdl    = Absorption model for WV (default 'ROS98')
%       absmdl.wvres = wv resonant absorption
%       absmdl.wvcnt = wv continuuum absorption
%     outputs: 
%          awet   = water vapor absorption profile (np/km)
%          adry   = dry air absorption profile (np/km)
%     subroutines:  
%          H2O_xxx = computes water vapor absorption
%          O2_xxx  = computes oxygen (dry air) absorption
%          (H2O and O2 both call subroutine Shape.)

function [awet,adry] = ClrAbs_xxx(p,tk,rho,frq,absmdl);

if nargin < 5
   absmdl.wvres = 'ROS98'; % default
   absmdl.wvcnt = 'ROS98'; % default
end

nl = length(p);
awet = zeros(size(p));
adry = zeros(size(p));
%factor = .182 * frq;
%db2np = log (10.) * 0.1;

for i = 1:nl
  
%c  Compute inverse temperature parameter; convert wet and dry p to kpa.
    %v = 300. / tk(i);
    %ekpa = e(i) / 10.;
    %pdrykpa = p(i) / 10. - ekpa;
  
%c  Compute H2O and O2 absorption (dB/km) and convert to np/km. 
    %[npp,ncpp] = H2O_xxx(pdrykpa,v,ekpa,frq);
    %awet(i) = (factor * (npp + ncpp)) * db2np;
    %[npp,ncpp] = O2N2_xxx(pdrykpa,v,ekpa,frq);
    %adry(i) = (factor * (npp + ncpp)) * db2np;
    % Per adesso uso le routine già adattate a partire dai .pro di Uli
    % in seguito al limite dovrò adattare quanto sopra
    ABSO2 = absO2R98(tk(i),p(i),rho(i),frq);
    ABSWV = absWVR98(tk(i),p(i),rho(i),frq,absmdl.wvres,absmdl.wvcnt);
    %ABSO2 = absO2_ARTS_PWR(tk(i),p(i),rho(i),frq,'98'); % NB!! Solo per prova...da togliere e rimettere le due linee sopra!!
    %ABSWV = absWV_ARTS_PWR(tk(i),p(i),rho(i),frq);      % NB!! Solo per prova...da togliere e rimettere le due linee sopra!!
    ABSN2 = absN2(tk(i),p(i),frq);
    awet(i) = ABSWV;
    adry(i) = ABSO2 + ABSN2;

end




return
%   EXPonential INTegration:  
%   Integrate the profile in array x over the layers defined in 
%   array ds, saving the integrals over each layer.
%     The algorithm assumes that x decays exponentially over each layer.
%
%     inputs passed as arguments: 
%	   zeroflg = flag to handle zero values (0:layer=0, 1:layer=avg)
%          x       = profile array
%          ds      = array of layer depths (km)
%          ibeg    = lower integration limit (profile level number)
%          iend    = upper integration limit (profile level number)
%	       factor  = factor by which result is multiplied (e.g., unit change)
%     outputs passed as arguments: 
%          xds     = array containing integrals over each layer ds
%          sxds    = integral of x*ds over levels ibeg to iend

function [sxds,xds] = ExpInt_xxx(zeroflg,x,ds,ibeg,iend,factor);

sxds = 0.0;
xds = zeros(size(ds));

for i = ibeg+1:iend
    
%c  Check for negative x value. If found, output message and return.
    if ((x(i-1) < 0.0) | (x(i) < 0.0))
        disp('Error encountered in ExpInt_xxx.m');
        return
%c  Find a layer value for x in cases where integration algorithm fails.
    elseif abs(x(i) - x(i-1)) < 1e-9
        xlayer = x(i);
    elseif (x(i-1) == 0.0) | (x(i) == 0.0)
        if zeroflg == 0
            xlayer = 0.0;
        else
            xlayer = (x(i) + x(i-1)) * 0.5;
        end
    else
%c      Find a layer value for x assuming exponential decay over the layer.
        xlayer = (x(i) - x(i-1)) / log(x(i) / x(i-1));
    end
        
%c  Integrate x over the layer and save the result in xds.
    xds(i) = xlayer * ds(i);
    sxds = sxds + xds(i);
        
end

sxds = sxds * factor;  
      
return
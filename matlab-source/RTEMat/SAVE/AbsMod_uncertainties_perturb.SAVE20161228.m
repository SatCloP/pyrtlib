% AMU = AbsMod_uncertainties_perturb(perturb_mode,perturb_what);
%
% perturb_mode: 'non' - no perturbation
%               'max' - original value + uncert           
%               'min' - original value - uncert           
%               'ran' - original value + uncert*random           
%
% perturb_what: 'all' or name of parameter as in AMU = AbsMod_uncertainties(0) 
%
% Example:
%          AMU = AbsMod_uncertainties_perturb('min','all');
%
% History:
% 2016/12 - Nico: first version

function AMU = AbsMod_uncertainties_perturb(perturb_mode,perturb_what);

% Set the spectroscopic parameters
dum = 0; % meaningless for the moment
AMU = AbsMod_uncertainties(dum);


% Select param to be perturbed
switch perturb_what
    case 'all'
        param = fieldnames(AMU);
        param = char(param);
    otherwise
        param = perturb_what;
end
npar = length(param(:,1));


% Perturb the spectroscopic parameters
for ipar = 1:npar
    
    param1 = deblank(param(ipar,:));
    p1 = getfield(AMU,param1); % get 1 parameter
    
    switch perturb_mode
        case 'non' % no perturbation            
            % Then don't do anything
            
        case 'max' % value + uncert           
            p1.value = p1.value + p1.uncer;
            
        case 'min' % value - uncert           
            p1.value = p1.value - p1.uncer;
            
        case 'ran' % value + uncert*random           
            % here rand or randn shall be used?
            % note that rand is only positive ([0-1]) and randn gives more
            % weight to values close to zero.
            % So, one way could be to combine them, using rand and the sign of randn
            % Assuming uncorrelated uncertainty:
            np1 = length(p1.value);
            for ip1 = 1:np1
                p1.value(ip1) = p1.value(ip1) + p1.uncer(ip1) * rand(1) * sign(randn(1));
            end
            
        otherwise % NaN to avoid misuse in case of typos           
            p1.value = NaN;
            
    end
    
    AMU = setfield(AMU,param1,p1); % set 1 parameter after perturbation

end

return 
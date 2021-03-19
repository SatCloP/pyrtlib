% AMU = AbsMod_uncertainties_perturb(perturb_mode,perturb_what,mdl);
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
% 2018/12 - Nico: modified to work for different models (e.g. ros98 has one additional line)

function AMU = AbsMod_uncertainties_perturb(perturb_mode,perturb_what,mdl);

if nargin < 3
   mdl = 0; % as it was for ACP 2018 paper
end

% Set the spectroscopic parameters
AMU = AbsMod_uncertainties(mdl);

% Select param to be perturbed
switch perturb_what
    case 'all'
        param = fieldnames(AMU);
        param = char(param);
    otherwise
        param = perturb_what;
        blank = strfind(param,' ');
        if ~isempty(blank)
           param_indx = str2num(param(blank+1:end));
           param = param(1:blank);
        else
           param_indx = []; 
        end
end
npar = length(param(:,1));


% Perturb the spectroscopic parameters
for ipar = 1:npar
    
    % get 1 parameter
    param1 = deblank(param(ipar,:));
    p1 = getfield(AMU,param1); 
    
    % in case param_indx is not set, take all indices
    if isempty(param_indx) 
        param_indx = 1:length(p1.value);
    end

    switch perturb_mode
        case 'non' % no perturbation            
            % Then don't do anything
            
        case 'max' % value + uncert           
            p1.value(param_indx) = p1.value(param_indx) + p1.uncer(param_indx);
            
        case 'min' % value - uncert           
            p1.value(param_indx) = p1.value(param_indx) - p1.uncer(param_indx);
            
        case 'ran' % value + uncert*random           
            % here rand or randn shall be used?
            % note that rand is only positive ([0-1]) and randn gives more
            % weight to values close to zero.
            % So, one way could be to combine them, using rand and the sign of randn
            % Assuming uncorrelated uncertainty:
            for ip1 = param_indx
                p1.value(ip1) = p1.value(ip1) + p1.uncer(ip1) * rand(1) * sign(randn(1));
            end
            
        otherwise % NaN to avoid misuse in case of typos           
            p1.value = NaN;
            
    end
    
    AMU = setfield(AMU,param1,p1); % set 1 parameter after perturbation

end

return 
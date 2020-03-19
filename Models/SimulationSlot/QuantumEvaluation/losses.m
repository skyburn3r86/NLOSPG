% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [field] = losses(field, z, varargin)
% Loss caller function. Bridging functionality and different models.
%   Arguments: Field, propagation distance z and type of propagation. Default: 'None'
    options = struct(...
                'type', 'None', ... 
                'zmax', 200e-6);
    optionNames=fieldnames(options);
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[])
        inpName = pair{1};
       if any(strcmp(inpName,optionNames))
           if isa(pair{2}, 'cell')
                options(1).(inpName) = pair{2}{1};
                options(2).(inpName) = pair{2}{2};
           else
               options(1).(inpName) = pair{2};
           end

       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

    if strcmp(options.('type'), 'None')
    elseif strcmp(options.('type'), 'Perturbation')
        field = PerturbationLosses(field, z/field.c_const);
    elseif strcmp(options.('type'), 'Scattering')
        field = ScatteringLosses(field, z);
    elseif strcmp(options.('type'), 'Exponential')
        field = ExponentialLosses(field, z);
    elseif strcmp(options.('type'), 'Beamsplitter')
        field = Beamsplitter(field, z, options.('zmax'));
    end

end
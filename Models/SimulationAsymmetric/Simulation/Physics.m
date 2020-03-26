% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Physics(model, Boundaries, Domains, varargin)
%PHYSICS Summary of this function goes here
%   Detailed explanation goes here
    options = struct(...
            'mat', 'Si', ...
            'V', {0, '[V]'});

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
               options(2).(inpName) = pair{2};
           end

       else
          error('%s is not a recognized parameter name',inpName)
       end
    end

    % Add ElectroDynamics Frequency Domain
    model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
%     model.component('comp1').physics('ewfd').selection.set(Domains.('halfDomain'));
    model.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);
    model.component('comp1').physics('ewfd').feature('sctr1').selection.set(Boundaries.('OuterBoundaries'));
%     model.component('comp1').physics('ewfd').create('PEC1', 'PerfectElectricConductor'); 
%     model.component('comp1').physics('ewfd').feature('PEC1').selection.set(Boundaries.('SymmetryLine'));
end

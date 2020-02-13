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
    model.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);
    model.component('comp1').physics('ewfd').feature('sctr1').selection.set(Boundaries.('OuterBoundaries'));

    % Add Graphene (Kubo Formulism)
    model.component('comp1').physics('ewfd').create('scu1', 'SurfaceCurrent', 1);
    model.component('comp1').physics('ewfd').feature('scu1').selection.set(Boundaries.('Graphene'));
    model.component('comp1').physics('ewfd').feature('scu1').set('Js0', {'sigmaxx*ewfd.Ex' 'sigmayy*ewfd.Ey' 'sigmazz*ewfd.Ez'});

    % Add Electrostatics
    model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
    model.component('comp1').physics('es').selection.set([Domains.('Organics') Domains.('Al2O3')]);
    % Create Ground Potential
    model.component('comp1').physics('es').create('gnd1', 'Ground', 1);
    model.component('comp1').physics('es').feature('gnd1').selection.set(Boundaries.('TopContact'));

    % Create Electrode Potential
    model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 1);
    model.component('comp1').physics('es').feature('pot1').selection.set(Boundaries.('Graphene'));
    % Set Voltage
    V = '0 [V]';
    if isa(options(1).('V'), 'double') || isa(options(1).('V'), 'float') || isa(options(1).('V'), 'int')
        V = [num2str(options(1).('V')) options(2).('V')];
    else
        V = options(1).('V');
    end
    model.component('comp1').physics('es').feature('pot1').set('V0', V);
end

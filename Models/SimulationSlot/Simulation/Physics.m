% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Physics(model, varargin)
%PHYSICS Summary of this function goes here
%   Detailed explanation goes here
    options = struct(...
            'mat', 'Si', ...
            'V', {0, '[V]'});
    ElectroStatics = false;
    optionNames=fieldnames(options);

    materials = varargin{1};
    materialNames=fieldnames(materials);

    % Add ElectroDynamics Frequency Domain
    model.component('comp1').physics.create('ewfd', 'ElectromagneticWavesFrequencyDomain', 'geom1');
    model.component('comp1').physics('ewfd').create('sctr1', 'Scattering', 1);
    % define here materials that touch the scattering boundary. Edges of this materail that
    % are not in contact with the simulation box are ignored.
    material_scattering = {'Substrate', 'OEO'};
    selection = [];
    for jj = 1:length(material_scattering)
        object = mphgetselection(model.selection(['geom1_' material_scattering{jj} '_bnd']));
        selection = [selection object.entities];
    end
    % selection cross-check via
    figure(1)
    mphgeom(model, 'geom1', 'Edgelabels', 'on');
    model.component('comp1').physics('ewfd').feature('sctr1').selection.set(selection);

    % Add Graphene (Kubo Formulism)
    model.component('comp1').physics('ewfd').create('scu1', 'SurfaceCurrent', 1);
    % define string label of graphene
    material_scattering = {'Graphene'};
    selection = [];
    for jj = 1:length(material_scattering)
        object = mphgetselection(model.selection(['geom1_' material_scattering{jj} '_bnd']));
        selection = [selection object.entities];
    end
    % selection cross-check via
    model.component('comp1').physics('ewfd').feature('scu1').selection.set(selection);
    model.component('comp1').physics('ewfd').feature('scu1').set('Js0', {'sigmaxx*ewfd.Ex' 'sigmayy*ewfd.Ey' 'sigmazz*ewfd.Ez'});

    % Add Electrostatics - here limited to the OEO domain between top and
    % bottom silicon
    if ElectroStatics
        model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
        material_scattering = {'OEO'};
        selection = [];
        for jj = 1:length(material_scattering)
            object = mphgetselection(model.selection(['geom1_' material_scattering{jj} '_dom']));
            selection = [selection object.entities];
        end
        model.component('comp1').physics('es').selection.set(selection);

        % Create Ground Potential
        model.component('comp1').physics('es').create('gnd1', 'Ground', 1);
        selection = [12 17];        % I do not like
        model.component('comp1').physics('es').feature('gnd1').selection.set(selection);

        % Create Electrode Potential
        model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 1);
        selection = [4 10];
        model.component('comp1').physics('es').feature('pot1').selection.set(selection);
        model.component('comp1').physics('es').feature('pot1').set('V0', 'V_bias');
    end
end

% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [model, Selections, BoundarySelection] = Geometry(model, varargin)
%GEOMETRY Generates Geometry of Initialized Parameters (Setup)
    % Get Arguments
    Selections = struct(...
        'SiO2', [], ...
        'SiNx', [], ...
        'Au', [], ...
        'Si', [], ...
        'Organics', [],...
        'Al2O3', [], ...
        'Air', []);

    global SimOps
    options = SimOps.getOptions;
    
    model.component('comp1').geom.create('geom1', 2);
    materialNames = fieldnames(Selections); 
    for jj = 1:length(materialNames)
        selection_container_dummy = model.component('comp1').geom('geom1').selection.create(materialNames{jj}, 'CumulativeSelection');
        selection_container_dummy.label(materialNames{jj});
    end

    % ***************** Create Geometry
        eps = 15e-9;
        counter = 1;

        % Substrate
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Thermal Oxide');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hSubstrate'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' '0'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'SiO2');
        counter = counter + 1;

        if isempty(Selections.('SiO2'))
            Selections.('SiO2') = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; -eps options(1).('hSubstrate')+eps], 'domain');
        else
            old = Selections.('SiO2');
            new = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; -eps options(1).('hSubstrate')+eps], 'domain');
            Selections.('SiO2') = [old, new];
        end
        if options(1).hRidge > 0
            % WG - base
            model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('BaseWaveguide');
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hRidge'});
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate'});
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', options(2).('mat'));
            counter = counter + 1;
        end

        % WG Right Ridge
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('WaveguideRRidge');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'(1-r)*wWG' 'hWG'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'wSim/2+(r-0.5)*wWG' 'hSubstrate+hRidge'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', options(2).('mat'));
        counter = counter + 1;

        % Organics, middle
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('OrganicsSlot');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'r*wWG' 'hWG'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'(wSim-wWG)/2' 'hSubstrate+hRidge'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Organics');
        counter = counter + 1;

        % Organics top
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('OrganicsTop');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hOrganic'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hRidge+hWG'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Organics');
        counter = counter + 1;

        % Cladding L
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('CladdingL');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'(wSim-wWG)/2' 'hWG'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hRidge'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'SiO2');
        counter = counter + 1;

        % Cladding R
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('CladdingR');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'(wSim-wWG)/2' 'hWG'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'(wSim+wWG)/2' 'hSubstrate+hRidge'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'SiO2');
        counter = counter + 1;

        % Air Buffer
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('AirTop');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hAir'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hRidge+hWG+hOrganic'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Air');
        counter = counter + 1;

        % Buffer
        if ~(options(1).('hBuffer') == 0)
            model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Buffer Oxide');
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hBuffer'});
            model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hOrganic'});
            counter = counter + 1;
        end

        outerBND = [];
        model.component('comp1').geom('geom1').run;
        for ii = 1:length(materialNames)
            temp = mphgetselection(model.selection(['geom1_' materialNames{ii} '_dom']));
            Selections.(materialNames{ii}) = temp.entities;
            temp = mphgetselection(model.selection(['geom1_' materialNames{ii} '_bnd']));
            outerBND = setxor(outerBND, temp.entities);
        end

        Graphene = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')-eps options(1).('hSubstrate')+eps], 'boundary');
        TopContact = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+eps], 'boundary');
        BoundarySelection = struct(...
            'OuterBoundaries', outerBND, ...
            'Graphene', Graphene,...
            'TopContact', TopContact);
end


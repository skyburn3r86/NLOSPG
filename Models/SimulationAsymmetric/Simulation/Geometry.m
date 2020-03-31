% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF

function [model, Selections, BoundarySelection] = Geometry(model, varargin)
%GEOMETRY Generates Geometry of Initialized Parameters (Setup)
    % Get Arguments
    materials = varargin{1};
    materialNames=fieldnames(materials);
    % ***************** Create Geometry
    model.component('comp1').geom.create('geom1', 2);

    % generting the selction containers for tracking materials of geometry
    % 1st entry lowest priority; last layer highest priority.
    % Higher priority: Material assignments of lower priority layers are overwritten
    for jj = 1:length(materialNames)
        selection_container_dummy = model.component('comp1').geom('geom1').selection.create(materialNames{jj}, 'CumulativeSelection');
        selection_container_dummy.label(materialNames(jj));
    end

    % ***************** Create Geometry
    eps = 15e-9;
    counter = 1;

    % Substrate
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('Thermal Oxide');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hSubstrate'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' '0'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Substrate');
    counter = counter + 1;

    % Parse parameter
    hRidge = char(model.param.get('hRidge'));
    if contains(hRidge, '[um]')
        hRidge = str2num(erase(hRidge, ' [um]'));
    elseif contains(hRidge, '[nm]')
        hRidge = str2num(erase(hRidge, ' [nm]'));
    elseif contains(hRidge, '[m]')
        hRidge = str2num(erase(hRidge, ' [m]'));
    else
        error('hRidge: Unit not recognized')
    end
    if hRidge > 0
        % WG - base
        model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('BaseWaveguide');
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hRidge'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate'});
        model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'PhotonicWG');
        counter = counter + 1;
    end

    % WG Right Ridge
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('WaveguideRRidge');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'(1-r)*wWG' 'hWG'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'wSim/2+(r-0.5)*wWG' 'hSubstrate+hRidge'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'PhotonicWG');
    counter = counter + 1;

    % Organics, middle
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('OrganicsSlot');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'r*wWG' 'hWG'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'(wSim-wWG)/2' 'hSubstrate+hRidge'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'OEO');
    counter = counter + 1;

    % Organics top
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('OrganicsTop');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hOrganic'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hRidge+hWG'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'OEO');
    counter = counter + 1;

    % Cladding L
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('CladdingL');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'(wSim-wWG)/2' 'hWG'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hRidge'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Substrate');
    counter = counter + 1;

    % Cladding R
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('CladdingR');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'(wSim-wWG)/2' 'hWG'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'(wSim+wWG)/2' 'hSubstrate+hRidge'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Substrate');
    counter = counter + 1;

    % Air Buffer
    model.component('comp1').geom('geom1').create(['r' num2str(counter)], 'Rectangle');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).label('AirTop');
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('size', {'wSim' 'hAir'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('pos', {'0' 'hSubstrate+hRidge+hWG+hOrganic'});
    model.component('comp1').geom('geom1').feature(['r' num2str(counter)]).set('contributeto', 'Air');
    counter = counter + 1;


    outerBND = [];
    % generating geometry
    model.component('comp1').geom('geom1').run;

    % Verifing the substrate selection and prints a report in the command line
    if(0)
        for jj = 1:length(materialNames)
            label2print = 'dom'; % Define here if you want to identifiy point 'pnt' or boundary 'bnd' or domain 'dom'
            result = mphgetselection(model.selection(['geom1_' materialNames{jj} '_' label2print]));
            materials.(materialNames{jj}) = result.entities;
            disp(['Layer' num2str(jj) '-' materialNames{jj} '_' label2print ': [' num2str(result.entities) ']']);
            pause(0.1);
        end
        disp('KEEP IN MIND: Higher layer number should overwrite lower layer number during Material ASSIGNMENT');
        % generates matlab plot of geometry and the domain, boundary, point
        % labels
        figure(1)
        switch label2print
            case 'dom'
                mphgeom(model, 'geom1', 'Facelabels', 'on');
            case 'bnd'
                mphgeom(model, 'geom1', 'Edgelabels', 'on');
            case 'pnt'
                mphgeom(model, 'geom1', 'vertexlabels', 'on');
            otherwise
        end
    end
end


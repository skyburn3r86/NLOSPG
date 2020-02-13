% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Geometry(model, varargin)
%GEOMETRY Generates Geometry of Initialized Parameters (Setup)
    % Get Arguments
    % Selections are containers for geometry entities of certain material
    
    materials = varargin{1}
    materialNames=fieldnames(materials);
    % ***************** Create Geometry
    model.component('comp1').geom.create('geom1', 2);
    
    % generting the selction containers for tracking materials of geometry
    % 1st entry lowest priority; last highest priority overriding lower
    % layers
    for jj = 1:length(materialNames)
        selection_container_dummy = model.component('comp1').geom('geom1').selection.create(materialNames{jj}, 'CumulativeSelection');
        selection_container_dummy.label(materialNames(jj));        
    end

    % Generating Objects: 
    
    % Generating substrate
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'ThermalOxide'], 'Rectangle');
    geom_dummy.label('Thermal Oxide');
	geom_dummy.set('base','center'); % base = center -> pos = center; base = corner -> pos = left bot corner of objects
    geom_dummy.set('pos', {'0' '-hSubstrate/2'}); % first position x, second position y
    geom_dummy.set('size', {'wSim' 'hSubstrate'});
    model.component('comp1').geom('geom1').feature('r_ThermalOxide').set('contributeto', materialNames(1));
    
    % generating cladding 
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'cladding'], 'Rectangle');
    geom_dummy.label('cladding');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' 'hcladding/2+hOEO'});
    geom_dummy.set('size', {'wSim' 'hcladding'});    
    model.component('comp1').geom('geom1').feature('r_cladding').set('contributeto', materialNames(1));
    
    % generating bottom WG burried in SiO2
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'photonic_wg_bot'], 'Rectangle');
    geom_dummy.label('photonic_wg_bot');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' '-hWG_bot/2'});
    geom_dummy.set('size', {'wWG' 'hWG_bot'});    
    model.component('comp1').geom('geom1').feature('r_photonic_wg_bot').set('contributeto', materialNames(2));

    % generating top WG 
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'photonic_wg_top'], 'Rectangle');
    geom_dummy.label('photonic_wg_top');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' 'hWG_top/2 + hOEO'});
    geom_dummy.set('size', {'wWG' 'hWG_top'});    
    model.component('comp1').geom('geom1').feature('r_photonic_wg_top').set('contributeto', materialNames(2));
           
    % generating OEO filling slot 
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'OEO_slot'], 'Rectangle');
    geom_dummy.label('OEO_slot');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' 'hOEO/2'});
    geom_dummy.set('size', {'wSim' 'hOEO'});    
    model.component('comp1').geom('geom1').feature('r_OEO_slot').set('contributeto', materialNames(3));
    
    % generating bottom graphene electrode entities
    geom_dummy = model.component('comp1').geom('geom1').create('poly_graphene_bot', 'Polygon');
    geom_dummy.set('source', 'vectors');
    % type: 'solid' filled object, 'closed' last and first point are conected, 'open' only lines between defined points of vector 
    geom_dummy.set('type', 'solid'); 
    geom_dummy.set('x', '-wSim/2 wWG/2');
    geom_dummy.set('y', '0 0');
    model.component('comp1').geom('geom1').feature('poly_graphene_bot').set('contributeto',  materialNames(6));
    
    % generating top graphene electrode entities
    geom_dummy = model.component('comp1').geom('geom1').create('poly_graphene_top', 'Polygon');
    geom_dummy.set('source', 'vectors');
    % type: 'solid' filled object, 'closed' last and first point are conected, 'open' only lines between defined points of vector 
    geom_dummy.set('type', 'solid'); 
    geom_dummy.set('x', '-wWG/2 wSim/2');
    geom_dummy.set('y', 'hOEO hOEO');
    model.component('comp1').geom('geom1').feature('poly_graphene_top').set('contributeto', materialNames(6));
    
    
    % generating geometry
    model.component('comp1').geom('geom1').run;
          
	% Verifing the substrate selection and prints a report in the command line  
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
    
    % alternative you can save the model and check the results
    % mphsave(model, 'test1234');
    
%     all = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; -eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+options(1).('hContact')+options(1).('hAir')+eps], 'boundary');
%     inner = mphselectbox(model, 'geom1', [+eps options(1).('wSim')-eps;eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+options(1).('hContact')+options(1).('hAir')-eps], 'boundary');
%     BoundarySelection = setdiff(all, inner);
%     Graphene = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')-eps options(1).('hSubstrate')+eps], 'boundary');
%     TopContact = mphselectbox(model, 'geom1', [-eps options(1).('wSim')+eps; options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')-eps options(1).('hSubstrate')+options(1).('hOrganic')+options(1).('hBuffer')+eps], 'boundary');
    BoundarySelection = []; % struct(...
%         'OuterBoundaries', BoundarySelection, ...
%         'Graphene', Graphene,...
%         'TopContact', TopContact);
end


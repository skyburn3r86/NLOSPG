% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function model = Geometry(model, varargin)
%GEOMETRY Generates Geometry of Initialized Parameters (Setup)
    % Get Arguments
    % Selections are containers for geometry entities of certain material
    
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

    % Generating Objects: 
    
    % Generating substrate
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'ThermalOxide'], 'Rectangle');
    geom_dummy.label('Thermal Oxide');
	geom_dummy.set('base','center'); % base = center -> pos = center; base = corner -> pos = left bot corner of objects
    geom_dummy.set('pos', {'0' '-hSubstrate/2'}); % first position x, second position y
    geom_dummy.set('size', {'wSim' 'hSubstrate'});    
    model.component('comp1').geom('geom1').feature('r_ThermalOxide').set('contributeto', materialNames{1});
    
    % generating cladding 
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'cladding'], 'Rectangle');
    geom_dummy.label('cladding');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' 'hcladding/2+hMetal/2+hOrganic/2'});
    geom_dummy.set('size', {'wSim' 'hcladding+hMetal+hOrganic'});    
    model.component('comp1').geom('geom1').feature('r_cladding').set('contributeto', materialNames(2));      
    
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'OEO_cladd'], 'Rectangle');
    geom_dummy.label('OEO_cladd_left');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' 'hOrganic/2'});
    geom_dummy.set('size', {'wSim' 'hOrganic'});    
    model.component('comp1').geom('geom1').feature('r_OEO_cladd').set('contributeto', materialNames(3));
              
    % generating OEO filling slot 
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'OEO_slot'], 'Rectangle');
    geom_dummy.label('OEO_slot');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'0' 'hMetal/2'});
    geom_dummy.set('size', {'wSlot' 'hMetal'});    
    model.component('comp1').geom('geom1').feature('r_OEO_slot').set('contributeto', materialNames(4));
    
%     % generating right Metal
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'metal_right'], 'Rectangle');
    geom_dummy.label('metal_right');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'-wMetal/2-wSlot/2' 'hMetal/2'});
    geom_dummy.set('size', {'wMetal' 'hMetal'});    
    model.component('comp1').geom('geom1').feature('r_metal_right').set('contributeto', materialNames(5));
    
%     % generating right Metal skindepth section
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'metal_right_skin'], 'Rectangle');
    geom_dummy.label('metal_right_skindepth');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'-skindepthMetal/2-wSlot/2' 'hMetal/2'});
    geom_dummy.set('size', {'skindepthMetal' 'hMetal'});    
    model.component('comp1').geom('geom1').feature('r_metal_right_skin').set('contributeto', materialNames(6));

    % generating left Metal 
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'metal_left'], 'Rectangle');
    geom_dummy.label('metal_left');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'wMetal/2+wSlot/2' 'hMetal/2'});
    geom_dummy.set('size', {'wMetal' 'hMetal'});    
    model.component('comp1').geom('geom1').feature('r_metal_left').set('contributeto', materialNames(5));
        
%     % generating right Metal skindepth section
    geom_dummy = model.component('comp1').geom('geom1').create(['r_' 'metal_left_skin'], 'Rectangle');
    geom_dummy.label('metal_left_skindepth');
	geom_dummy.set('base','center');
    geom_dummy.set('pos', {'skindepthMetal/2+wSlot/2' 'hMetal/2'});
    geom_dummy.set('size', {'skindepthMetal' 'hMetal'});    
    model.component('comp1').geom('geom1').feature('r_metal_left_skin').set('contributeto', materialNames(6));
    
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


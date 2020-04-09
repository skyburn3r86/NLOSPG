% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = meshing(model, varargin)
%MESH Meshes the structure. Can be given arguments: Materials 
%'mesh[mat]' for the meshing parameters of the material, 'mesh[mat]x' and
%'mesh[mat]y' for anisotropic meshing (x, y=1 will be isotropic)
    % Define Material Parameters
    materials = varargin{1};
    materialNames =fieldnames(materials);

    % A good starting point for a mesh shize depend on the character of the
    % mode. A plasmonic mesh should be finder then a photonic to resolve
    % the exponential decay in the metal (Screening length ie 50nm for 1550nm)
    % A rule of thumb for photonic mesh size is wavelength (wl) over 8
    % times the refractive index for photonics, while plasmonics should be
    % on the order of < ScreeningLength/5. In doubt perform a mesh analysis study by
    % sweeping the mesh size and tracking the solution of the effective
    % refractive index to see convergence of the solution for reduced mesh
    % size within the comsol model saved by mphsave(model, 'Convergence
    % Test')
    
    model.component('comp1').mesh.create('mesh1');
    selected_domains = [];
    % reverse sweep direction to give high priority of meshing for later material
    for jj = length(materialNames):-1:1
        % checks if material properties are defined by txt file
       
        if ~isempty(strfind(materials.(materialNames{jj}), '.txt'))
            % check if plasmonic waveguide 
            if ~isempty(strfind(lower(materialNames{jj}), 'metalfine'))
                meshsize = 'wl/plasmonic_mesh'; % plasmonic_mesh is defined in ModelSetup_Parameters
                % reduce simulation time by utilzing scaling factors. For
                % instance, hybrid waveguide with metal extending towards
                % inifity for y --> yscale < 1 increases meshsize;
                xscale = 1;
                yscale = 1;
            elseif ~isempty(strfind(lower(materialNames{jj}), 'electrodes')) || ~isempty(strfind(lower(materialNames{jj}), 'metalrough'))
                meshsize = 'wl/20';                
                xscale = 1;
                yscale = 1;
            elseif ~isempty(strfind(lower(materialNames{jj}), 'substrate')) || ~isempty(strfind(lower(materialNames{jj}), 'cladding'))
                meshsize = 'wl/10';                
                xscale = 1;
                yscale = 1;
            elseif ~isempty(strfind(lower(materialNames{jj}), 'wg'))
                refractive_index = real(extractRefractiveIndex(materials.(materialNames{jj}), 'model', model));
                meshsize = ['wl/' num2str(refractive_index) '/8/3'];               
                xscale = 1;
                yscale = 1;
            else % case of photonic waveguide/components
                % open the txt file of the material to extract the data and
                % interpolate the refractive index
                refractive_index = real(extractRefractiveIndex(materials.(materialNames{jj}), 'model', model));
                meshsize = ['wl/' num2str(refractive_index) '/8/3'];       
                xscale = 1;
                yscale = 1;
            end
            % adding triangular mesh
            model.component('comp1').mesh('mesh1').create(['ftri', materialNames{jj}], 'FreeTri');
            model.component('comp1').mesh('mesh1').feature(['ftri', materialNames{jj}]).label(['ftri', materialNames{jj}]);
            model.component('comp1').mesh('mesh1').feature(['ftri', materialNames{jj}]).set('xscale', xscale);
            model.component('comp1').mesh('mesh1').feature(['ftri', materialNames{jj}]).set('yscale', yscale);
            model.component('comp1').mesh('mesh1').feature(['ftri', materialNames{jj}]).selection.geom('geom1', 2);
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).create('size1', 'Size');
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).feature('size1').set('hauto', 1);
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).feature('size1').set('custom', 'on');
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).feature('size1').set('hmax', meshsize);
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).feature('size1').set('hmaxactive', true);
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).feature('size1').set('hmin', 3.86E-11);
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).feature('size1').set('hminactive', false);
            
            % adding domains to mesh materialNames{jj}
            objects = mphgetselection(model.selection(['geom1_' materialNames{jj} '_dom']));
            mesh_selection = [];
            % for-loop scans if the current domain has been asigned
            % previously. 
            for ii = 1:length(objects.entities)
                if isempty(find(selected_domains == objects.entities(ii)))
                    mesh_selection = [mesh_selection objects.entities(ii)];
                    selected_domains = [selected_domains objects.entities(ii)];
                end
            end
            model.component('comp1').mesh('mesh1').feature(['ftri',  materialNames{jj}]).selection.set(mesh_selection);
        end  
    end
    model.component('comp1').mesh('mesh1').run;
end


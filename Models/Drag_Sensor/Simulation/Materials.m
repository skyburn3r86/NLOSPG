% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Materials(model,varargin)
%MATERIALS Creates the Materials for the Model.

    % Define Material Parameters
    materials = varargin{1};
    materialNames=fieldnames(materials);
	global library_path 
    
    % Define material properties for ElectroDynamics argument
    for jj = 1:length(materialNames)
        %         if ~isempty(strfind(materials.(materialNames{jj}), '.txt'))
        material_forloop = materials.(materialNames{jj});
        switch material_forloop(1:3)
            case 'nk_' % optical properties
                models.(['n_' materialNames{jj}]) = ['n' materialNames{jj} '(wl*1e9)'];
                models.(['k_' materialNames{jj}]) = ['k' materialNames{jj} '(wl*1e9)'];
                material_forloop = strrep(material_forloop,'nk_', 'eps_');
            case 'ana'
                % TO DO
            otherwise
                models.(['n_' materialNames{jj}]) = '0';
                models.(['k_' materialNames{jj}]) = '1';
        end
    end
    
    for ii = 1:length(materialNames)
        % set optical properties. 
        if ~isempty(strfind(materials.(materialNames{ii}), '.txt'))
            n = models.(['n_' materialNames{ii}]);
            k = models.(['k_' materialNames{ii}]);
            model_dummy = model.component('comp1').material.create(['mat' materialNames{ii}], 'Common');
            model_dummy.propertyGroup.create('RefractiveIndex', 'Refractive index');
            model_dummy.label(materialNames{ii});
            model_dummy.propertyGroup('RefractiveIndex').set('n', {n '0' '0' '0' n' '0' '0' '0' n});
            model_dummy.propertyGroup('RefractiveIndex').set('ki', {k '0' '0' '0' k '0' '0' '0' k});
            objects = mphgetselection(model.selection(['geom1_' materialNames{ii} '_dom']));
            model_dummy.selection.set(objects.entities);            
            model.component('comp1').material(['mat' materialNames{ii}]).propertyGroup('def').set('relpermittivity', {['eps_' materialNames{ii}]});
            
            % checks if rf files have been initalized
            material_forloop = materials.(materialNames{ii});
            if 2 == exist([library_path '\DataIn\' material_forloop]) % 0 does not exist, 2 file exist
                eps_str = ['eps' materialNames{ii} '_re(f_rf)+i*eps' materialNames{ii} '_im(f_rf)'];
                model.component('comp1').material(['mat' materialNames{ii}]).propertyGroup('def').set('relpermittivity', {eps_str});
            end
            
            % Adapt here material properties by the user given values
            if ~isempty(strfind(lower(materialNames{ii}),'high_k')) % || ~isempty(strfind(lower(materialNames{ii}),'oeo'))
                try
                    % checks if parameter is defined. If not error message.
                    model.param.get(['eps_' materialNames{ii}]);
                    model.component('comp1').material(['mat' materialNames{ii}]).propertyGroup('def').set('relpermittivity', {['eps_' materialNames{ii}]});
                catch
                    error_prompt = ['Error in ComsolEvaluation/Materials parameter eps_' materialNames{ii} ' is not defined. Add Parameter or remove if case'];
                    error(error_prompt);
                end
            end
        end
    end
end



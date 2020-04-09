% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [model] = Materials(model,varargin)
%MATERIALS Creates the Materials for the Model.

    % Define Material Parameters
    materials = varargin{1};
    materialNames=fieldnames(materials);
    
    % Define material properties for ElectroDynamics argument 
    for jj = 1:length(materialNames)
        if ~isempty(strfind(materials.(materialNames{jj}), '.txt'))
            if strcmp(materialNames{jj}, 'OEOWG')
                % special case template for utilizing birefringent metals.
                models.(['n_' materialNames{jj} '_eoAxis']) = ['real(sqrt(eps' materialNames{jj} '_re(wl*1e9)+i*eps' materialNames{jj} '_im(wl*1e9)))'];
                models.(['k_' materialNames{jj} '_eoAxis']) = ['imag(sqrt(eps' materialNames{jj} '_re(wl*1e9)+i*eps' materialNames{jj} '_im(wl*1e9)))'];
            end
            models.(['n_' materialNames{jj}]) = ['real(sqrt(eps' materialNames{jj} '_re(wl*1e9)+i*eps' materialNames{jj} '_im(wl*1e9)))'];
            models.(['k_' materialNames{jj}]) = ['imag(sqrt(eps' materialNames{jj} '_re(wl*1e9)+i*eps' materialNames{jj} '_im(wl*1e9)))'];
        % model specific part
        elseif ~isempty(strfind(materials.(materialNames{jj}), 'High_k'))
            models.(['n_' materialNames{jj}]) = 'n_High_k';
            models.(['k_' materialNames{jj}]) = '0';            
        else
            models.(['n_' materialNames{jj}]) = '1';
            models.(['k_' materialNames{jj}]) = '0';
        end
        
        % @TODO: if case for n and k files or eps.... 
    end

    for ii = 1:length(materialNames)
        name = materialNames{ii};
        if strcmp(name, 'OEOWG')
            no = models.(['n_' name]);
            ne = models.(['n_' name '_eoAxis']);
        else            
            no = models.(['n_' name]);
            ne = no;
        end
            k = models.(['k_' name]);
        model_dummy = model.component('comp1').material.create(['mat' name], 'Common');
        model_dummy.propertyGroup.create('RefractiveIndex', 'Refractive index');
        model_dummy.label(name);
%         model_dummy.propertyGroup('RefractiveIndex').set('n', '');
%         model_dummy.propertyGroup('RefractiveIndex').set('ki', '');
        model_dummy.propertyGroup('RefractiveIndex').set('n', {ne '0' '0' '0' no' '0' '0' '0' no});
        model_dummy.propertyGroup('RefractiveIndex').set('ki', {k '0' '0' '0' k '0' '0' '0' k});
        objects = mphgetselection(model.selection(['geom1_' materialNames{ii} '_dom']));
        model_dummy.selection.set(objects.entities);
           
        if ~isempty(strfind(lower(name),'high_k')) || ~isempty(strfind(lower(name),'oeo'))
            model.component('comp1').material(['mat' name]).propertyGroup('def').set('relpermittivity', {['eps_' name]});
        end
    end
end


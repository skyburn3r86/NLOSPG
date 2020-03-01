modelpath = pwd;
% path for global variabls, paths, material data, etc
cd('../');
initPaths(modelpath);

% Note, defaults value are set in ModelSetup_Parameters
% Sweep parameters. String has to match parameter name of the comsol model 
para_sweep{1}.values = linspace(10, 200, 21)*1e-9;
para_sweep{1}.str = 'hOEO';
para_sweep{1}.unit = '[m]';
para_sweep{2}.values = linspace(100, 250, 21)*1e-9;
para_sweep{2}.str = 'hWG_bot';
para_sweep{2}.unit = '[m]';
para_sweep{3}.values = linspace(1550, 1900, 3)*1e-9;
para_sweep{3}.str = 'wl';
para_sweep{3}.unit = '[m]';

% maping N-dimensional parameter sweep onto linear list
[param_list] = combParameterSweep(para_sweep);

for idx_param_list = 1:size(param_list.values,1)
    % displays current run
    disp(['Current run:' num2str(round(idx_param_list/size(param_list.values,1)*100,100)) '%  --- '...
        param_list.print{idx_param_list}]);    
    % ModelSetup_Parameters - init Model and defines parameters    
    [comsol_model, materials, sim_parameters] = ModelSetup_Parameters(param_list, idx_param_list,...
        'n_start', {3.4, ' '});
    % setup model
    comsol_model = Geometry(comsol_model, materials);
    comsol_model = Materials(comsol_model, materials);
    comsol_model = meshing(comsol_model, materials);
    comsol_model = Physics(comsol_model, materials);
    try
        % calculate compute
        [comsol_model] = Compute(comsol_model);
    catch
        disp(['Error at Simulating' param_list.print{idx_param_list}])
    end
    
    % Save & Evaluate
    if 1
        save_str = strrep(param_list.print{idx_param_list},'>_<', '');
        save_str = strrep(save_str,'_','');
        save_str = strrep(save_str,'[','');
        save_str = strrep(save_str,']','');
        save_str = strrep(save_str,' ','_');
        mphsave(comsol_model, ['./ComsolModels/' save_str '.mph']);
    end
    sim_results{idx_param_list,1} = comsolEvaluation(comsol_model, sim_parameters, materials, 'title', save_str);
end
%%
if 1
    % to do save as h5 file
    workspace_save_str = '__SIS__';
    for dimension = 1:length(para_sweep)
        workspace_save_str = [workspace_save_str strrep(para_sweep{dimension}.str,'_','') '_'];
    end
    save(['./Results/' workspace_save_str '_Results.mat']);
end

plotResultList(sim_results, param_list,'data_str', 'g_0',...
    'para_str', {para_sweep{1}.str, 'hWG_bot' ,'wl'}, 'para_values', {[], [] ,1550e-9},'data_FOM', 'max');

%% Data Evaluation
for idx_row = 1:length(hOEO_array)
    for idx_col = 1:length(hWG_bot_array)
        g_0_plot(idx_row,idx_col) = sim_results{idx_row,idx_col}.g_0.value;
        gamma_plot(idx_row,idx_col) = sim_results{idx_row,idx_col}.gamma.value;
        Q_plot(idx_row,idx_col) = sim_results{idx_row,idx_col}.Q.value;    
        neff_plot(idx_row,idx_col) = sim_results{idx_row,idx_col}.neff.value; 
        ng_plot(idx_row,idx_col) = sim_results{idx_row,idx_col}.ng.value;              
    end
end

close all
plotResultList(result_List, para_list, '')  

figure(5)
surface(hWG_bot_array, hOEO_array, real(neff_plot))
colorbar
xlabel('hWG bot [m]')
ylabel('h OEO [m]')
figure(6)
surface(hWG_bot_array, hOEO_array, Q_plot )
colorbar
xlabel('hWG bot [m]')
ylabel('h OEO [m]')


figure(7)
title('n_{group} [m]')
surface(hWG_bot_array, hOEO_array, ng_plot)
colorbar
xlabel('hWG bot [m]')
ylabel('h OEO [m]')

figure(8)
title('Vacuum Cooperetivity')
surface(hWG_bot_array, hOEO_array, real(g_0_plot).^2./1e6./193e12.*Q_plot)
colorbar
xlabel('hWG bot [m]')
ylabel('h OEO [m]')
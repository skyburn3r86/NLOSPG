modelpath = pwd;
% path for global variabls, material data, etc
initPaths('Models/SimulationSlot');

para_sweep{1}.values = linspace(10, 200, 21)*1e-9;
para_sweep{1}.str = 'hOEO';
para_sweep{1}.unit = '[m]';
para_sweep{2}.values = linspace(50, 250, 9)*1e-9;
para_sweep{2}.str = 'dSlot';
para_sweep{2}.unit = '[m]';
para_sweep{3}.values = linspace(1550, 1900, 2)*1e-9;
para_sweep{3}.str = 'wl';
para_sweep{3}.unit = '[m]';
para_sweep{3}.values = linspace(300, 900, 3)*1e-9;
para_sweep{3}.str = 'wWG';
para_sweep{3}.unit = '[m]';

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
        continue
    end
    
    % Save & Evaluate
    save_str = strrep(param_list.print{idx_param_list},'>_<', '');
    save_str = strrep(save_str,'_','');
    save_str = strrep(save_str,'[','');
    save_str = strrep(save_str,']','');
    save_str = strrep(save_str,' ','_');
    if 1
        mphsave(comsol_model, ['./ComsolModels/' save_str '.mph']);
    end
    sim_results{idx_param_list,1} = comsolEvaluation(comsol_model, sim_parameters, materials, 'title', save_str);
end
if 1
    % to do save as h5 file
    workspace_save_str = '__SIS__';
    for dimension = 1:length(para_sweep)
        workspace_save_str = [workspace_save_str strrep(para_sweep{dimension}.str,'_','') '_'];
    end
    save(['./Results/' workspace_save_str '_Results.mat']);
end

%% Data Evaluation
close all
% availabe simulatin results
available_sim_results = [];
for jj = 1:length(sim_results{1})
    available_sim_results = [available_sim_results '  /  ' sim_results{1}(jj).str];
end
display(['Available simulation results:' available_sim_results]);

% availabe simulatin results
available_param = [];
for jj = 1:length(param_list.str)
    available_param = [available_param '  /  ' param_list.str{jj}];
end
display(['Parameters:' available_param]);

%%
plotResultList(sim_results, param_list,'data_str', 'g_0',...
    'para_str', {para_sweep{1}.str, 'hWG_bot' ,'wl'}, 'para_values', {[], [] ,1550e-9},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);

plotResultList(sim_results, param_list,'data_str', 'g_0',...
    'para_str', {para_sweep{1}.str, 'hWG_bot' ,'wl'}, 'para_values', {[], [] ,1900e-9},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);
%%
plotResultList(sim_results, param_list,'data_str', 'gamma',...
    'para_str', {para_sweep{1}.str, 'hWG_bot' ,'wl'}, 'para_values', {[], [] ,1550e-9},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);

plotResultList(sim_results, param_list,'data_str', 'gamma',...
    'para_str', {para_sweep{1}.str, 'hWG_bot' ,'wl'}, 'para_values', {[], [] ,1900e-9},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);

%%
plotResultList(sim_results, param_list,'data_str', 'Q', ...
    'para_str', {para_sweep{1}.str, 'hWG_bot' ,'wl'}, 'para_values', {[], [] ,1550e-9},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);

%%
plotResultList(sim_results, param_list,'data_str', 'Q', ...
    'para_str', {'wl', para_sweep{1}.str}, 'para_values', {[], []},...
    'data_FOM', 'max', 'save', ['./Results/' workspace_save_str]);
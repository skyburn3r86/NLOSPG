clear all
close all

modelpath = pwd;
% path for global variabls, paths, material data, etc
cd('../');
cd('../');
initPaths(modelpath);

% Note, defaults value are set in ModelSetup_Parameters
% Sweep parameters. String has to match parameter name of the comsol model 
para_sweep{1}.values = linspace(30, 300, 8)*1e-9;
para_sweep{1}.str = 'wOEO';
para_sweep{1}.unit = '[m]';
para_sweep{2}.values = linspace(50, 200, 4)*1e-9;
para_sweep{2}.str = 'hOEO';
para_sweep{2}.unit = '[m]';
para_sweep{3}.values = [5 50 500];
para_sweep{3}.str = 'eps_High_k';
para_sweep{3}.unit = '';
para_sweep{4}.values = linspace(100, 300, 3)*1e-9;
para_sweep{4}.str = 'hWG';
para_sweep{4}.unit = '[m]';
para_sweep{5}.values = linspace(400, 800, 3)*1e-9;
para_sweep{5}.str = 'wWG';
para_sweep{5}.unit = '[m]';
para_sweep{6}.values = linspace(2000, 3000, 2)*1e-9;
para_sweep{6}.str = 'wElectrode';
para_sweep{6}.unit = '[m]';


% maping N-dimensional parameter sweep onto linear list
[param_list] = combParameterSweep(para_sweep);

% defining save location
save_folder = generateSavePath(para_sweep);

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
    
    % Save Comsol Model & Evaluate
        file_str = strrep(param_list.print{idx_param_list},'>_<', '');
        file_str = strrep(file_str,'_','');
        file_str = strrep(file_str,'[','');
        file_str = strrep(file_str,']','');
        file_str = strrep(file_str,' ','_');
    if 1
        try
            mphsave(comsol_model, [save_folder '/ComsolModels/' file_str '.mph']);
        catch
        end
    end
    sim_results{idx_param_list,1} = comsolEvaluation(comsol_model, sim_parameters, materials, 'title', file_str, 'path', save_folder);   
end

%% Saving results as jason files -> move to seperate function!

save([save_folder '\rawData.mat']);
% adapt to also save sim_parameters and param_sweep
            writeToJson(param_list, sim_results, 'results_and_paramters')



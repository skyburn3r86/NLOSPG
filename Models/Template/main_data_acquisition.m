% clear all
close all

modelpath = pwd;
% path for global variabls, paths, material data, etc
cd('../');
cd('../');
initPaths(modelpath);

% Note, defaults value are set in ModelSetup_Parameters
% Sweep parameters. String has to match parameter name of the comsol model 
para_sweep{1}.values = linspace(10, 200, 21)*1e-9;
para_sweep{1}.str = 'hOEO';
para_sweep{1}.unit = '[m]';
para_sweep{2}.values = linspace(100, 250, 7)*1e-9;
para_sweep{2}.str = 'hWG_bot';
para_sweep{2}.unit = '[m]';
para_sweep{3}.values = linspace(1550, 1900, 1)*1e-9;
para_sweep{3}.str = 'wl';
para_sweep{3}.unit = '[m]';
para_sweep{4}.values = linspace(300, 900, 3)*1e-9;
para_sweep{4}.str = 'wWG';
para_sweep{4}.unit = '[m]';

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
    
    % Save Comsol Model & Evaluate
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

%% Saving results as jason files -> move to seperate function!
    % generating folder name
    save_str = [];
    for jj = 1:length(para_sweep)
        if length(para_sweep{jj}.values) > 1
            save_str = [save_str '__' para_sweep{jj}.str ...
                '_(' num2str(min(para_sweep{jj}.values)) '-' num2str(max(para_sweep{jj}.values)) ')'...
                para_sweep{jj}.unit];
        else
            save_str = [save_str '__' para_sweep{jj}.str '_(' num2str(min(para_sweep{jj}.values)) ')'...
                para_sweep{jj}.unit];
        end
    end
    save_str = strrep(save_str , '.', 'pt');
    % generating folder
    [status, msg, msgID] = mkdir(['./Results\' save_str])
    % check if folder already exist. if true ask user if override should be
    % performed
    if strcmp(msg,'Directory already exists.')
        prompt = {'Folder already exist. Options: Redfine path string, overwrite (if name is not changed) or cancel saving:'};
        title = 'File Path';
        dims = [1 length(save_str)+20];
        definput = {save_str};
        save_str = inputdlg(prompt,title,dims,definput); % cancel will return empty str
        if ~isempty(save_str)
            [status, msg, msgID] = mkdir(['./Results\' save_str])
            save(['./Results/' workspace_save_str '_Results.mat']);
            % adapt to also save sim_parameters and param_sweep
%             writeToJson(param_list, sim_results, 'results_and_paramters')           
        end
    else
        save(['./Results/' workspace_save_str '_Results.mat']);
        % adapt to also save sim_parameters and param_sweep
%         writeToJson(param_list, sim_results, 'results_and_paramters')
    end
 


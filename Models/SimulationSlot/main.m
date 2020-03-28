clear all; 
modelpath = pwd;
% path for global variabls, material data, etc
cd('../..');
initPaths(modelpath);


para_sweep{1}.values = linspace(100, 200, 3)*1e-9;
para_sweep{1}.str = 'hWG';
para_sweep{1}.unit = '[m]';
para_sweep{2}.values = linspace(25, 125, 5)*1e-9; 
para_sweep{2}.str = 'dSlot';
para_sweep{2}.unit = '[m]';
para_sweep{3}.values = linspace(450, 800, (800-450)/25 + 1)*1e-9;
para_sweep{3}.str = 'wWG';
para_sweep{3}.unit = '[m]';

global old_Ep
global old_neff

old_neff = 2.18; 
old_Ep = 0; 
hWG = 220; 
[param_list] = combParameterSweep(para_sweep);

for idx_param_list = 1:size(param_list.values,1)
    % displays current run
    disp(['Current run:' num2str(round(idx_param_list/size(param_list.values,1)*100,100)) '%  --- '...
        param_list.print{idx_param_list}]);    

    % ModelSetup_Parameters - init Model and defines parameters    
    [comsol_model, materials, sim_parameters] = ModelSetup_Parameters(param_list, idx_param_list,...
        'n_start', {old_neff, ' '}, 'hWG', {hWG, '[nm]'});
    
    % setup model
    comsol_model = Geometry(comsol_model, materials);
    comsol_model = Materials(comsol_model, materials);
    comsol_model = meshing(comsol_model, materials);
    comsol_model = Physics(comsol_model, materials);
    
    try
        % calculate compute
        [comsol_model] = Compute(comsol_model, 'mat', 'Si');
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
    if 0
        mphsave(comsol_model, ['./ComsolModels/' save_str '.mph']);
    end
    sim_results{idx_param_list,1} = comsolEvaluation(comsol_model, sim_parameters, materials, 'title', save_str);
end
JsonName = ['./Results/Data__hWG-' num2str(hWG) 'nm_Mode20.json'];
writeToJson(param_list, sim_results, JsonName); 


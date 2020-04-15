clear all; 
modelpath = pwd;
% path for global variabls, material data, etc
cd('../..');
initPaths(modelpath);


para_sweep{3}.values = linspace(100, 200, 2)*1e-9;
para_sweep{3}.str = 'hOEO';
para_sweep{3}.unit = '[m]';
para_sweep{2}.values = linspace(75, 125, 3)*1e-9; 
para_sweep{2}.str = 'dSlot';
para_sweep{2}.unit = '[m]';
para_sweep{1}.values = linspace(700, 1100, (1100-600)/50 + 1)*1e-9;
para_sweep{1}.str = 'wWG';
para_sweep{1}.unit = '[m]';

global old_Ep
global old_neff

old_neff = 1.5671; 
old_Ep = 0; 
hWG = 800; 

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
    if 1
        mphsave(comsol_model, ['./ComsolModels/' save_str '.mph']);
    end
    [sim_results{idx_param_list,1}, mode_results{idx_param_list, 1}] = comsolEvaluation(comsol_model, sim_parameters, materials, 'title', save_str, 'expression', 'real(-i*ewfd.Ex)');
end
JsonName = ['./Results/Data__hWG-' num2str(hWG) 'nm_Mode10.json'];
writeToJson(param_list, sim_results, JsonName, 'model', comsol_model); 
JsonName = ['./Results/Modes__hWG-' num2str(hWG) 'nm_Mode10.json'];
writeToJson(param_list, mode_results, JsonName, 'model', comsol_model, 'sollabels', {'wl_1310_', 'wl_2620_'}); 

modelpath = pwd;
% path for global variabls, material data, etc
initPaths('Models/SimulationSlot');

% para_sweep{1}.values = linspace(150, 200, 1)*1e-9;
para_sweep{1}.values = 150*1e-9;
para_sweep{1}.str = 'hOEO';
para_sweep{1}.unit = '[m]';
para_sweep{2}.values = 100e-9; 
para_sweep{2}.str = 'dSlot';
para_sweep{2}.unit = '[m]';
para_sweep{3}.values = linspace(700, 800, 101)*1e-9;
para_sweep{3}.str = 'wWG';
para_sweep{3}.unit = '[m]';
hWGA = [340]; 
modes = [2.72];
modeOrders = [2]; 
global old_Ep
global old_neff
global order
for idx_hWG = 1:length(hWGA)
    for idx_Mode = 1:size(modes, 2)
        hWG = hWGA(idx_hWG); 
        mode = modes(idx_hWG, idx_Mode); 
        old_neff = mode; 
        old_Ep = 0; 
        order = modeOrders(idx_hWG, idx_Mode); 
        [param_list] = combParameterSweep(para_sweep);
        clear simresults;
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
        JsonName = ['./Results/Data__hWG-' num2str(hWGA(idx_hWG)) 'nm_Mode-' num2str(idx_Mode) '.json'];
        writeToJson(param_list, sim_results, JsonName); 
    end
end

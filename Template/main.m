modelpath = pwd;
% path for global variabls, paths, material data, etc
cd('../');
initPaths(modelpath);

% define all variables in standard units using e-X to avoid confusion of
% units. Is it nm is it microns...?
wWGA = 600*1e-9;
% loop1_str
hOEO_array = linspace(10, 200, 11)*1e-9;
hWG_bot_array =  linspace(10, 200, 11)*1e-9;
wl = 1550e-9;
% reads refractive index of photonic waveguide as start value
n_start = real(extractRefractiveIndex(materials.PhotonicWG, 'wavelength', wl)); 
nr_modes = 4;

nu = cell(length(hOEO_array), length(hWG_bot_array));
losses = cell(length(hOEO_array), length(hWG_bot_array));

for idx_row = 1:length(hOEO_array)
    for idx_col = 1:length(hWG_bot_array)
        hWG_bot = hWG_bot_array(idx_col);
        hOEO = hOEO_array(idx_row);
        DumpingName = ['Row ' num2str(idx_row) '_Col ' num2str(idx_col)...
            '_hOEO ' num2str(hOEO) 'm_ hWG' num2str(hWG_bot) 'm_'];
        % Simulate
        disp(['Current Simulation:'...
            ' - Column_' num2str(round(idx_col/length(hWG_bot_array)*100)) '%' ...
            ' - Row_' num2str(round(idx_row/length(hOEO_array)*100)) '%']);
        % ModelSetup_Parameters - init Model and defines parameters
        [comsol_model, materials, sim_parameters] = ModelSetup_Parameters('wl', {wl, '[m]'}, 'n_start', {n_start, ' '},...
            'hOEO', {hOEO, '[m]'}, 'hWG_bot', {hWG_bot, '[m]'}, 'nr_modes', {nr_modes, ' '});
        sim_parameters(1).idx_row = idx_row;        
        sim_parameters(1).idx_col = idx_col;
        % TODO DispersionRelation - what is the purpose of this? Shouldnt we
        % include that into our material libary? 
%         [Simulation, refr] = DispersionRelation(Simulation, 'plot', false); 
        % Builiding the geometry based on ModelSetup_Parameters
        comsol_model = Geometry(comsol_model, materials);
        comsol_model = Materials(comsol_model, materials); 
        comsol_model = meshing(comsol_model, materials); 
        comsol_model = Physics(comsol_model, materials); 
        try
            [comsol_model] = Compute(comsol_model); 
        catch
            disp(['Error at Simulating' DumpingName])
        end

        % Save & Evaluate
        if 1
            mphsave(comsol_model, ['./ComsolModels/' DumpingName '_Model.mph']);
        end
        sim_results{idx_row, idx_col} = comsolEvaluation(comsol_model, sim_parameters, materials);
%         save(['./ComsolModels/' DumpingName '_Results.mat'], 'Results');
    end
end

PlotLosses(losses)
PlotEfficiciency(hOEO,OrganicA, nu)
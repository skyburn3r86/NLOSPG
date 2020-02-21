modelpath = pwd;
% path for global variabls, paths, material data, etc
cd('../');
initPaths(modelpath);

% define all variables in standard units using e-X to avoid confusion of
% units. Is it nm is it microns...?
wWGA = 600*1e-9;
hWGA_array = [220, 260, 340, 400, 500]*1e-9;
hOrganicA = linspace(200, 300, 3)*1e-9; 
wl = 1550e-9;
% reads refractive index of photonic waveguide as start value
n_start = real(extractRefractiveIndex(materials.PhotonicWG, 'wavelength', wl)); 
nr_modes = 4;

nu = cell(length(hWGA_array), length(hOrganicA));
losses = cell(length(hWGA_array), length(hOrganicA));

for idx_row = 1:length(hWGA_array)
    for idx_col = 1:length(hOrganicA)
        hOrganic = hOrganicA(idx_col);
        hWG = hWGA_array(idx_row);
        DumpingName = [num2str(hOrganic) 'm_' num2str(hWG) 'm_'];
        % Simulate
        % ModelSetup_Parameters - init Model and defines parameters
        [comsol_model, materials, sim_parameters] = ModelSetup_Parameters('wl', {wl, '[m]'}, 'n_start', {n_start, ' '},...
            'hWG_top', {hWG, '[m]'}, 'hOrganic', {hOrganic, '[m]'}, 'hBuffer', {0, '[m]'},...
        'nr_modes', {nr_modes, ' '});
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
        if 0
            mphsave(comsol_model, ['./ComsolModels/' DumpingName '_Model.mph']);
        end
        [Results] = comsolEvaluation(comsol_model, sim_parameters, materials);
        save(['./ComsolModels/' DumpingName '_Results.mat'], 'Results');
        %dumpStruct(Results, ['./Data/200109_' DumpingName '_ResultsDump.csv']);
%         try
%             PlotDispersionRelation(Results, 'fileName', DumpingName);
%         catch
%            disp(['Error at Plotting' DumpingName])
%            continue
%         end

        nu{idx_row, idx_col} = QuantumEvaluation(comsol_model, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        losses{idx_row, idx_col} = AssignResults(Results, hWG, hOrganic);
    end
end

PlotLosses(losses)
PlotEfficiciency(hWGA_array,OrganicA, nu)
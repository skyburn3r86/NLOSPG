modelpath = pwd;
% path for global variabls, material data, etc
global library_path 
cd('../');
addpath([library_path '\DataIn'])
cd('DataIn');
library_path = pwd;
cd(modelpath)
addpath([modelpath '\Evaluation'])
addpath([modelpath '\Plotting'])
addpath([modelpath '\QuantumEvaluation'])
addpath([modelpath '\Simulation'])

N = 25;
wWGA = 600;
hWGA_array = [220, 260, 340, 400, 500];
hOrganicA = linspace(200, 300, 3); 
material = 'Si';            %Si or SiNx

nu = cell(length(hWGA_array), length(hOrganicA));
losses = cell(length(hWGA_array), length(hOrganicA));

for ii = 1:length(hWGA_array)
    for kk = 1:length(hOrganicA)
        hOrganic = hOrganicA(kk);
        hWG = hWGA_array(ii);
        DumpingName = [num2str(hOrganic) 'nm_' num2str(hWG) 'nm_'];
        % Simulate
        % ModelSetup_Parameters - init Model and defines parameters
        [comsol_model, materials] = ModelSetup_Parameters('Nl', {N, ' '}, 'hWG_top', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        % DispersionRelation - what is the purpose of this? Shouldnt we
        % include that into our material libary? 
%         [Simulation, refr] = DispersionRelation(Simulation, 'plot', false); 
        % Builiding the geometry based on ModelSetup_Parameters
        comsol_model = Geometry(comsol_model, materials);
        
%         [Simulation, Selections, BoundarySel] = Geometry(Simulation, 'mat', {0, material}, 'Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hBuffer', {0, '[nm]'});
        [comsol_model] = Materials(comsol_model, materials); 
        [comsol_model] = meshing(comsol_model, Selections, 'mat', material); 
        [comsol_model] = Physics(comsol_model, BoundarySel, Selections); 
        try
            [comsol_model] = Compute(comsol_model, 'mat', material); 
        catch
            disp(['Error at Simulating' DumpingName])
        end

        % Save & Evaluate
        mphsave(comsol_model, ['./ComsolModels/' DumpingName '_Model.mph']);
        [Results] = Evaluation(comsol_model, 'Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        save(['./ComsolModels/' DumpingName '_Results.mat'], 'Results');
        %dumpStruct(Results, ['./Data/200109_' DumpingName '_ResultsDump.csv']);
        try
            PlotDispersionRelation(Results, 'fileName', DumpingName);
        catch
           disp(['Error at Plotting' DumpingName])
           continue
        end

        nu{ii, kk} = QuantumEvaluation(comsol_model, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        losses{ii, kk} = AssignResults(Results, hWG, hOrganic);
    end
end

PlotLosses(losses)
PlotEfficiciency(hWGA_array,OrganicA, nu)
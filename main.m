N = 25;
wWGA = 600;
hWGA = [220, 260, 340, 400, 500];
hOrganicA = linspace(200, 300, 3); 
material = 'Si';            %Si or SiNx

nu = cell(length(hWGA), length(hOrganicA));
losses = cell(length(hWGA), length(hOrganicA));

for ii = 1:length(hWGA)
    for kk = 1:length(hOrganicA)
        hOrganic = hOrganicA(kk);
        hWG = hWGA(ii);
        DumpingName = [num2str(hOrganic) 'nm_' num2str(hWG) 'nm_'];
        % Simulate
        Simulation2q = Setup('Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        [Simulation, refr] = DispersionRelation(Simulation, 'plot', false); 
        [Simulation, Selections, BoundarySel] = Geometry(Simulation, 'mat', {0, material}, 'Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hBuffer', {0, '[nm]'});
        [Simulation] = Materials(Simulation, Selections); 
        [Simulation] = meshing(Simulation, Selections, 'mat', material); 
        [Simulation] = Physics(Simulation, BoundarySel, Selections); 
        try
            [Simulation] = Compute(Simulation, 'mat', material); 
        catch
            disp(['Error at Simulating' DumpingName])
        end

        % Save & Evaluate
        mphsave(Simulation, ['./ComsolModels/' DumpingName '_Model.mph']);
        [Results] = Evaluation(Simulation, 'Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        save(['./ComsolModels/' DumpingName '_Results.mat'], 'Results');
        %dumpStruct(Results, ['./Data/200109_' DumpingName '_ResultsDump.csv']);
        try
            PlotDispersionRelation(Results, 'fileName', DumpingName);
        catch
           disp(['Error at Plotting' DumpingName])
           continue
        end

        nu{ii, kk} = QuantumEvaluation(Simulation, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hOrganic', {hOrganic, '[nm]'}, 'hBuffer', {0, '[nm]'});
        losses{ii, kk} = AssignResults(Results, hWG, hOrganic);
    end
end

PlotLosses(losses)
PlotEfficiciency(hWGA,OrganicA, nu)
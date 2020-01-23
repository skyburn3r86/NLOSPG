N = 25;
wWGA = 400; 
hWGA = [220, 340];
material = 'Si';            %Si or SiNx
QSStates = cell(length(wWGA), 1); 
for ii = 1:length(wWGA)
    for kk = 1:length(hWGA)
        wWG = wWGA(ii);
        hWG = hWGA(kk); 
        DumpingName = [num2str(wWG) 'nm_' num2str(hWG) 'nm_'];
        % Simulate
        Simulation = Setup('Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hBuffer', {0, '[nm]'});
        [Simulation, refr] = DispersionRelation(Simulation, 'plot', false); 
        [Simulation, Selections, BoundarySel] = Geometry(Simulation, 'mat', {0, material}, 'Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'}, 'hBuffer', {0, '[nm]'});
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
        [Results] = Evaluation(Simulation, 'Nl', {N, ' '}, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'});
        save(['./ComsolModels/' DumpingName '_Results.mat'], 'Results');
        dumpStruct(Results, ['./Data/200109_' DumpingName '_ResultsDump.csv']);
        try
            PlotDispersionRelation(Results, 'fileName', DumpingName);
        catch
           disp(['Error at Plotting' DumpingName])
           continue
        end
        mainTest(Simulation, 'wWG', {wWG, '[nm]'}, 'hWG', {hWG, '[nm]'});
        disp(['Done with Quantum Evaluation ' num2str(ii/length(wWGA)*100) '%'])
    end
end
save(['./ComsolModels/' DumpingName '_QS.mat'], 'QSStates');

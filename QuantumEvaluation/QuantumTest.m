% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
%% Read Model, Read Ep, Es, Ei
model = mphload('C:\Users\kkeller\Documents\PhD\04 Projects\SinglePhotonGeneration\03 Simulation\01 Comsol\200123_NLOOpt\ComsolModels\200nm_220nm__Model.mph');

% Initialize the propagation distance and vector
zmin = 0;   
zmax = 50e-6; 
dz = 100e-9;
zz= linspace(0, zmax, (zmax-zmin)/dz - 1); 

[Ep, Es, Ei] = test_ExtractField(model);

% Read the first Solution, i.e. Nondegenerate
Ep = Ep{1};
Es = Es{1}; 
Ei = Ei{1}; 

% Normalize the Fields
Ep = Ep.normalizeField(zmax);
Es = Es.normalizeField(zmax); 
Ei = Ei.normalizeField(zmax); 

% Initialize Quantum States
Ep.NPhotons = 100; 
Es.NPhotons = 10; 
Ei.NPhotons = 10; 
Es = Es.initializeState(1); 
Ei = Ei.initializeState(1); 
Ep = Ep.initializeState(Ep.NPhotons + 1); 

%%
close all;

record = true;

if record
    writerObj = VideoWriter('Videos/Generated_withoutLoss.avi');
    writerObj.FrameRate = 15;
    open(writerObj);
    myfig = figure(1);
    hold on; 
end
[P, Estemp, Eitemp, Eptemp] =  QuantumCalc(Es, Ep, Ei, zz);
disp('Done')
for ii = 1:length(zz)
    if record
        plot(linspace(0, Estemp{ii}.NPhotons, Estemp{ii}.NPhotons + 1), conj(Estemp{ii}.psi).*Estemp{ii}.psi)
        tex = sprintf('z=%.1f', zz(ii)*1e6);
        text(0.8*Es.NPhotons, 0.9, tex) 
        ylim([0, 1]);
        xlabel('Photon Number')
        ylabel('Probability')
        drawnow
        myframe = getframe(gca);
%         size(myframe.cdata); 
        writeVideo(writerObj, myframe);
        clf(myfig);
    end
end

if record
    hold off;
    close(writerObj);
end

% Plot the power
figure(1)
semilogy(zz, real(P)/max(P))
ylim([1e-2, 1.1])

disp('Done')

% fileID = fopen(['Data/exp_' num2str(N(kk)) '.txt'], 'w');
% for ii = 1:length(z)
%     x = z(ii);
%     y = P(ii)/max(P);
%     fprintf(fileID, '%.7f\t%.14f\n', x, y);
% end
% fclose(fileID);


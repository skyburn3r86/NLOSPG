% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
%% Read Model, Read Ep, Es, Ei
model = mphload('C:\Users\kkeller\Documents\PhD\04 Projects\SinglePhotonGeneration\03 Simulation\01 Comsol\200123_NLOOpt\ComsolModels\200nm_220nm__Model.mph');

% Initialize the propagation distance and vector
zmin = 0;   
zmax = 200e-6; 
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
Es.NPhotons = 10; 
Ei.NPhotons = 10; 
Es = Es.initializeState(1); 
Ei = Ei.initializeState(1); 

% Initialize Pump beam
tau = 1e-6;                 % Measure during 1 microsecond
Pin = 1;                    % 1 W input Power
Ep = Ep.setTime(tau);       % Set duration of measurement
Ep = Ep.setPower(Pin);      % Set Input power

%% 
lossStr = 'Lossless';
if strcmp(lossStr, 'Lossless')
    chi1is = 0; 
    chi1ii = 0;
elseif strcmp(lossStr, 'Lossy')
    chi1is = -0.0172; 
    chi1ii = -0.0220;
else
    chi1is = Es.chi1i; 
    chi1ii = Ei.chi1i;
end
Es.chi1i = chi1is;
Ei.chi1i = chi1ii; 
%% Simulate
[P, Estemp, Eitemp, Eptemp] =  QuantumCalc(Es, Ep, Ei, zz);
disp('Done with Simulation')

%% Plot
close all;
record = true;
if record
    writerObj = VideoWriter('Videos/Generated_withoutLoss.avi');
    writerObj.FrameRate = 30;
    open(writerObj);
    myfig = figure(1);
    hold on; 
end
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

disp('Done Plotting')

%%

fileID = fopen(['Data/exp.txt'], 'w');
for ii = 1:length(zz)
    x = zz(ii);
    y = P(ii)/max(P);
    fprintf(fileID, '%.7f\t%.14f\n', x, y);
end
fclose(fileID);


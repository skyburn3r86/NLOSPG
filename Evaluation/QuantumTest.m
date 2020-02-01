% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
%% Read Model, Read Ep, Es, Ei
model = mphload('C:\Users\kkeller\Documents\PhD\04 Projects\SinglePhotonGeneration\03 Simulation\01 Comsol\200123_NLOOpt\ComsolModels\200nm_220nm__Model.mph');
% Coupled WaveEquations
%   Detailed Summary Here
z = linspace(0, 0.1e-3, 1000);
P = zeros(length(z), 1);
[Ep, Es, Ei] = test_ExtractField(model);
Ep = Ep{1};
Ep = Ep.normalizeField;
Es = Es{1}; 
Ei = Ei{1}; 
Es = Es.initializeState(1); 
Ei = Ei.initializeState(1); 
Ep.NPhotons = 100; 
Ep = Ep.initializeState(Ep.NPhotons + 1); 
%%
Ep.NPhotons = 100; 
Ep = Ep.initializeState(Ep.NPhotons + 1); 
zmin = 0;   
zmax = 100e-6; 
dz = 30e-6/51;
zz= linspace(0, zmax, (zmax-zmin)/dz - 1); 
P = zeros(length(zz), 1); 
% Create Movie from the electrical field

numberOfFrames = length(zz);

t = linspace(0, 5, numberOfFrames);
hFigure = figure;

% Set up the movie structure.
% Preallocate movie, which will be an array of structures.
% First get a cell array with all the frames.
vidHeight = 344;
vidWidth = 446;
record = true;
writerObj = VideoWriter('myVideo.avi');
writerObj.FrameRate = 15;
myfig = figure(1);
hold on; 
open(writerObj);
for ii = 1:length(zz)
    [g0, Pp, Ex, Ey, Ez, Hx, Hy, Hz] =  QuantumCalc(Es, Ep, Ei, zz(1:ii));
    P(ii) = Pp; 
    if record
        surf(real(Hz))
        zlim([-1.5e4, 1.5e4]);
        drawnow
        myframe = getframe(gca);
        size(myframe.cdata); 
        writeVideo(writerObj, myframe);
        clf(myfig);
    end
end
hold off;
z = zz; 
close(writerObj);

figure(1)
semilogy(zz, real(P)/max(P))
ylim([1e-2, 1.1])

fileID = fopen('exp.txt', 'w');
for ii = 1:length(z)
    x = z(ii);
    y = P(ii)/max(P);
    fprintf(fileID, '%.7f\t%.14f\n', x, y);
end
fclose(fileID);
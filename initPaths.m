% Generated through Matlab
% Author:           Christian Haffner 
% E-Mail:           christian.haffner1986@gmail.com
% Organization:     ETHZ ITET IEF/NIST
% Date 2020/02/14


function initPaths(modelpath)
%initPaths - loads the path of the base model, functions, class and material data
    global library_path;
    library_path = pwd;
    %% add here cross user files, class, 
    addpath([library_path '\DataIn']);
    addpath([library_path '\Classes']);
    addpath([library_path '\HelperFunctions']);
    
%     cd(modelpath)                 % Removed because it changes directory
%     to the modelpath and therefore the prefix modelpath in the following
%     sections throws an error!
    
    %% add here model specific folders
    addpath([modelpath '\Evaluation']);
    addpath([modelpath '\Evaluation']);
    addpath([modelpath '\Plotting']);
    addpath([modelpath '\QuantumEvaluation']);
    addpath([modelpath '\Simulation']); 
end

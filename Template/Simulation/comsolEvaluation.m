% Generated through Matlab
% Author:           Killian Keller
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function simulation_results = comsolEvaluation(model, simulation_parameters, materials, varargin)
%EVALUATION Evaluates the Simulation and can be adapted to your individual
%style
% VARARGIN options (extend on need)
% 1. 'expression' - expression of the comsol results to plott and save
% 2. 'nr_solution' - flag to save or  of the comsol results to plott and save
% 3. 'title' - string of the figure title will be expanded by neff
% 3. 'path' - string of the path figure will be stored to 

    warning('off', 'all')
    
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('Arguments needs propertyName/propertyValue pairs')
    end
    
    % finding modes
    [neffTE, nr_solutionTE, neffTM, nr_solutionTM] = findGuidedModes(model,...
        'substrate', materials.Substrate,'deltaN_threshold', 0.05, 'polarization_threshold', 0.5);
    
    % plotting and saving modes
    for jj = 1:length(nr_solutionTM)
        if ~isempty(nr_solutionTM)
            title_str = ['Row ' num2str(simulation_parameters(1).idx_row) ...
                '__Col ' num2str(simulation_parameters(1).idx_col) ...
                '__TMpol'];
            saveSolutionSnapshot(model, 'expression', 'ewfd.normE', 'nr_solution', nr_solutionTM(jj),...
                'title', title_str);
        end
    end
    
    % Calculations on the desired mode(s) - in this example it is the fundamental TM Mode
    simulation_results = calculateVaccumCoupling(model, 'active_material', 'OEO', 'nr_solution', nr_solutionTM(1));
         
        
% % % % % % % % %         Quantum Calculations..

% Normalization of the fields... 
        
end


% Generated through Matlab
% Author/last edit:           Killian Keller/Christian 
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function simulation_results = comsolEvaluation(model, simulation_parameters, materials, varargin)
%EVALUATION Evaluates the Simulation and can be adapted to your individual
%style 
%OUPUT: simulaiton_results - structure value feautring 1D array, str name of
%property

%Input: simulaiton_results
% VARARGIN options (extend on need)
% 1. 'expression' - expression of the comsol results to plott and save
% 2. 'title' - string of the figure title will be expanded by neff
% 3. 'path' - string of the path figure will be stored to 

    warning('off', 'all')
    
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('Arguments needs propertyName/propertyValue pairs')
    end
    expression_str = 'ewfd.normE'; 
    % transfroming variable inputs into variables
    for ii = 1:length(varargin)-1
        switch varargin{ii}
            case 'title'
                title_str = varargin{ii+1};
            case 'expression'
                expression_str = varargin{ii + 1}; 
            case 'path'
        end
    end
        
    % finding modes
    [neffTE, nr_solutionTE, neffTM, nr_solutionTM] = findGuidedModes(model,...
        'deltaN_threshold', 0.03, 'polarization_threshold', 0.5, 'OuterSolNum', [1, 2]);
    
    if isa(neffTE, 'cell')
        OuterSolNum = [1, 2]; 
        for ii = 1:length(nr_solutionTE)
            nr_solutionTEtemp = nr_solutionTE{ii}; 
            for jj = 1:length(nr_solutionTEtemp)
                if ~isempty(nr_solutionTEtemp)
                    saveSolutionSnapshot(model, 'expression', expression_str, 'nr_solution', nr_solutionTEtemp(jj),...
                        'title', [num2str(OuterSolNum(ii)) title_str], 'OuterSolNum', OuterSolNum(ii), 'dset', 'dset2');
                end
            end
        end
    else
        for jj = 1:length(nr_solutionTE)
            if ~isempty(nr_solutionTE)
                saveSolutionSnapshot(model, 'expression', expression_str, 'nr_solution', nr_solutionTE(jj),...
                    'title', title_str);
            end
        end
    end
    % plotting and saving modes
    global old_neff
    [~, I] = min(abs(neffTE{1} - old_neff)); 
    old_neff = neffTE{1}(I); 
    
    % Calculations on the desired mode(s) - in this example it is the fundamental TM Mode
   simulation_results = calculateVaccumCoupling(model, 'active_material', 'OEO', ...
       'nr_solution', [nr_solutionTE{1}(I), nr_solutionTE{2}(1)],'OuterSolNums', OuterSolNum, 'type', 'SPDC');
         
        
% % % % % % % % %         Quantum Calculations..

% Normalization of the fields... 
        
end


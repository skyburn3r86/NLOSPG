% Generated through Matlab
% Author/last edit:           Killian Keller/Christian 
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [simulation_results, error_prompt]= comsolEvaluation(model, simulation_parameters, materials, varargin)
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
    % default values
    error_prompt = [];
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
        error('Arguments needs propertyName/propertyValue pairs')
    end
    
    % transfroming variable inputs into variables
    for ii = 1:2:length(varargin)-1
        switch varargin{ii}
            case 'title'
                title_str = varargin{ii+1};
            case 'expression'
            case 'path'
                save_folder = varargin{ii+1};
        end
    end
        
    % finding modes
    [neffTE, nr_solutionTE, neffTM, nr_solutionTM] = findGuidedModes(model,...
        'substrate', materials.Substrate,'deltaN_threshold', 0.01, 'polarization_threshold', 0.5);
    if isempty(nr_solutionTM)
        error_prompt = 'Error in ComsolEvaluation/findGuidedModes no TM mode found';
        % Calculations on the desired mode(s) - in this example it is the fundamental TM Mode
        simulation_results = calculateVaccumCoupling(model, 'active_material', 'OEOWG', 'nr_solution', 1, 'error_flag', 1);
    else        
        % plotting and saving modes
        for jj = 1:length(nr_solutionTM(1))
            if ~isempty(nr_solutionTM)
                saveSolutionSnapshot(model, 'expression', 'abs(ewfd.Ez)', 'nr_solution', nr_solutionTM(jj),...
                    'title', title_str, 'path', save_folder);
            end
        end
        % Calculations on the desired mode(s) - in this example it is the fundamental TM Mode
        simulation_results = calculateVaccumCoupling(model, 'active_material', 'OEOWG', 'nr_solution', nr_solutionTM(1),...
            'expr_field1', 'ewfd.Ey', 'expr_field2', 'ewfd.Ey', 'expr_field3', 'es.Ey');        
        if isnan(simulation_results(1).value)
            error_prompt = 'Error in ComsolEvaluation/calculateVaccumCoupling while extracting solutions';
        end  
    end
    if ~isempty(error_prompt)
            errorProtocol(save_folder, strrep(title_str,'pt','.'), error_prompt);
    end
end



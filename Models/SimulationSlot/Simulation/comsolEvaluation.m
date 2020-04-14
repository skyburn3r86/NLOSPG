% Generated through Matlab
% Author/last edit:           Killian Keller/Christian 
% E-Mail:           killian.keller@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [simulation_results, mode_results] = comsolEvaluation(model, simulation_parameters, materials, varargin)
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
    
    % if no prior field has been calculated, do comparison with neff
    global old_Ep        
    global old_neff
    
    % Checking if the Simulation is the first one. If it is, compare with
    % the effective refractive index and read the intial field. If it is
    % not, compare the overlaps with the existing fields.
    if ~isa(old_Ep, 'ClassicalField')
        [~, I] = min(abs(neffTE{1} - old_neff)); 
        [~, old_Ep] = overlapFields(model, old_Ep, 'dset', 'dset2', 'N', 200, 'OuterSolNums', 1, 'SolNums', nr_solutionTE{1}(I));
    else
        overlap = zeros(1, length(nr_solutionTE{1})); 
        field = cell(1, length(nr_solutionTE{1})); 
        for idx_mode_TE = 1:length(nr_solutionTE{1})
           [overlap(idx_mode_TE), field{idx_mode_TE}] = overlapFields(model, old_Ep, 'dset', 'dset2', 'OuterSolNums', 1, 'N', 200, 'SolNums', nr_solutionTE{1}(idx_mode_TE)); 
        end
       [m, I] = max(overlap); 
        if m < 0.5
           % TODO: Implement trying further to correct the overlap
            disp('Warning:    overlap smaller than 50%! Keeping old mode as reference');
            mphsave(model); 
            neffStr = ''; 
            for idx_mode = 1:length(nr_solutionTE{1})
               neffStr = [neffStr num2str(idx_mode) ': ' num2str(real(neffTE{1}(idx_mode)), 4) ';  '];
            end
            prompt = {['Overlap not clear. Please select the correct mode! ' neffStr]};
            dims = [1 length(prompt{1})];
            title = 'Sim Index';
            definput = {'1'};
            I = inputdlg(prompt,title,dims,definput); % cancel will return empty str
            I = str2num(I{1}); 
            old_Ep = field{I};  
        else
           old_Ep = field{I};  
       end
       old_Ep = field{I};        
    end
    old_neff = neffTE{1}(I); 
    
    mode_results = struct();
    mode_results(1).str = 'TE';
    mode_results(1).value = neffTE; 
    mode_results(1).unit = '[a.u.]';
    mode_results(2).str = 'TM'; 
    mode_results(2).value = neffTM; 
    mode_results(2).unit = '[a.u.]';
    
    % Calculations on the desired mode(s) - in this example it is the fundamental TM Mode
   simulation_results = calculateVaccumCoupling(model, 'active_material', 'OEO', ...
       'nr_solution', [nr_solutionTE{1}(I), nr_solutionTE{2}(1)],'OuterSolNums', OuterSolNum, 'type', 'SPDC');
         
        
% % % % % % % % %         Quantum Calculations..

% Normalization of the fields... 
        
end


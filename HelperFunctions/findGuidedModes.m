% Generated through Matlab
% Author:           Christian Haffner
% E-Mail:           christian.haffner@ief.ee.ethz.ch
% Organization:     ETHZ ITET IEF
    
function [neffTE, nr_solutionTE, neffTM, nr_solutionTM] = findGuidedModes(model, varargin)
%Input Variables: - model = comsol model
% VARARGIN options (extend on need)
% 1. 'substrate' - name of the materials.Substrate txt file. This file
% needs to be defined in the DataIn folder of the framework
% 2. 'polarization_threshold' - minimal ratio between tranverse E-Field to determine TE or TM polarization
% 3. 'deltaN_threshold' - numerical value representing guided mode cutoff
% 4. 'energy_confinement_thresh' - numerical value representing desired 

% TO DO - add the node finder to check for order of mode
%       - power confinement criteria

    global library_path
    options = struct('OuterSolNum', 1); 
    %% checks that paris of input argements are defined
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('Arguments needs propertyName/propertyValue pairs')
    end

    %% default value definition
    % Threshold value (0.1) needs to be tested for you case.
    deltaN_threshold = 0.1;
    % Note, neff < ncladding/substrate -> no guided/physical mode.
    n_cladding = 1.44;
    polarization_threshold = 0.5;
    
    %% scans through varagin. -1 and +1 of for loop due to option/value pairs
    for ii = 1:length(varargin)-1
        switch varargin{ii}
            case 'substrate'
                n_cladding = real(extractRefractiveIndex(varargin{ii+1}, 'model', model)); 
            case 'deltaN_threshold'
                deltaN_threshold = varargin{ii+1};
            case 'polarization_threshold'
                polarization_threshold = varargin{ii+1};
            case 'OuterSolNum'
                options.OuterSolNum = varargin{ii+1};
            otherwise
        end
    end  
    dset = 'dset1'; 
    if length(options.OuterSolNum) > 1
        dset = 'dset2'; 
    end
    neffTE = cell(length(options.OuterSolNum), 1); 
    neffTM = cell(length(options.OuterSolNum), 1); 
    nr_solutionTM = cell(length(options.OuterSolNum), 1); 
    nr_solutionTE = cell(length(options.OuterSolNum), 1); 
    
    for ii = 1:length(options.OuterSolNum)
        %% % Evaluate Data, Sorting and Mode Selection
        % Find physical solutions
        % 1. Get the effective refractive indeces
        % dsetX - Dataset Tag read by comosl model -> right click solX -> Properties -> Tag
        % X is defined by the various studies. Here, X = 1 is electromagnetic (ewfd) mode solution, X = 2 is
        % 2 (es). 'outersolnum' is defined by the prameter sweep.
        neff = mphglobal(model, 'ewfd.neff', 'dataset', dset, 'outersolnum', options.OuterSolNum(ii));

        % sorting refractive index largest (fundamental) to smallest (higher order modes) value.
        idxPhysicalMode = (neff > n_cladding + deltaN_threshold);

        % Next differentiation between TE/TM polarization (measured parallel to substrate) is still needed.
        fieldEx = mphint2(model, 'abs(ewfd.Ex)^2', 'surface', 'dataset', dset, 'outersolnum', options.OuterSolNum(ii));
        fieldEy = mphint2(model, 'abs(ewfd.Ey)^2', 'surface', 'dataset', dset, 'outersolnum', options.OuterSolNum(ii));
        lumericalTETMdefinition = fieldEx./(fieldEy+fieldEx);

        % ratio > 2 indicates TE polarization - no guarantee. The threshold
        % value is an number based on experience and should not be take for
        % granted. Larger values -> stronger polarized mode.
        % TE Mode
        idxTE = (lumericalTETMdefinition > polarization_threshold);
        idxTE = idxTE' & idxPhysicalMode;
        % sorting the solutions following mode order - 1st = Fundamental
        counter = 1;
        nr_solutionTEtemp = [];
        neffTEtemp =[];
        for jj = 1:length(idxTE)
            if real(idxTE(jj)) > 0
                neffTEtemp(counter) = neff(jj);
                nr_solutionTEtemp = [nr_solutionTEtemp jj];
                counter = counter+1;
            end
        end
        [~, index] = sort(neffTEtemp,'descend');
        neffTE{ii} = neffTEtemp(index);
        nr_solutionTE{ii} = nr_solutionTEtemp(index);
        
        % TM Mode
        idxTM = (lumericalTETMdefinition < polarization_threshold);
        idxTM = idxTM' & idxPhysicalMode;
        counter = 1;
        nr_solutionTMtemp = [];
        neffTMtemp =[];
        for jj = 1:length(idxTM)
            if real(idxTM(jj)) > 0
                neffTMtemp(counter) = neff(jj);
                nr_solutionTMtemp = [nr_solutionTMtemp jj];
                counter = counter+1;
            end
        end
        [~, index] = sort(neffTMtemp,'descend');
        % Write Solutions
        neffTM{ii} = neffTMtemp(index);
        nr_solutionTM{ii} = nr_solutionTMtemp(index);

        % Alternatively, modes can be found by calculating the energy
        % confinement to the photonic waveguide
    end
    
    if length(options.OuterSolNum) == 1
       neffTE = neffTE{1}; 
       neffTM = neffTM{1}; 
       nr_solutionTM = nr_solutionTM{1};
       nr_solutionTE = nr_solutionTE{1}; 
    end
       
end